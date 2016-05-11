#include "fabm_driver.h"
#include "fabm_private.h"

module test_models

use fabm_types

implicit none

integer, parameter :: interior_state_offset = 0
integer, parameter :: surface_state_offset  = 1000
integer, parameter :: bottom_state_offset   = 2000

integer, parameter :: interior_dependency_offset = 3000
integer, parameter :: horizontal_dependency_offset = 4000

type,extends(type_base_model) :: type_test_model
   type (type_state_variable_id),        allocatable :: id_state(:)
   type (type_surface_state_variable_id),allocatable :: id_surface_state(:)
   type (type_bottom_state_variable_id), allocatable :: id_bottom_state(:)

   type (type_dependency_id)             :: id_dep
   type (type_dependency_id)             :: id_depth
   type (type_horizontal_dependency_id)  :: id_hz_dep

   integer :: nstate         = 12
   integer :: nsurface_state = 11
   integer :: nbottom_state  = 10
contains
   procedure :: initialize
   procedure :: do_surface
   procedure :: do_bottom
end type

   contains

subroutine initialize(self,configunit)
   class (type_test_model), intent(inout), target :: self
   integer,                 intent(in)            :: configunit

   integer          :: i
   character(len=8) :: strindex

   allocate(self%id_state(self%nstate))
   do i=1,self%nstate
      write (strindex,'(i0)') i
      call self%register_state_variable(self%id_state(i),'state'//trim(strindex),'','state variable #'//trim(strindex))
   end do
   allocate(self%id_surface_state(self%nsurface_state))
   do i=1,self%nsurface_state
      write (strindex,'(i0)') i
      call self%register_state_variable(self%id_surface_state(i),'surface_state'//trim(strindex),'','surface state variable #'//trim(strindex))
   end do
   allocate(self%id_bottom_state(self%nbottom_state))
   do i=1,self%nbottom_state
      write (strindex,'(i0)') i
      call self%register_state_variable(self%id_bottom_state(i),'bottom_state'//trim(strindex),'','bottom state variable #'//trim(strindex))
   end do
   call self%register_dependency(self%id_dep,standard_variables%temperature)
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_hz_dep,standard_variables%wind_speed)
end subroutine initialize

subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_test_model),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   integer  :: i
   real(rk) :: value

   _HORIZONTAL_LOOP_BEGIN_
      do i=1,self%nstate
         _GET_(self%id_state(i),value)
         if (value/=i+interior_state_offset) call self%fatal_error('do_surface','invalid value of interior state variable.')
         _SET_SURFACE_EXCHANGE_(self%id_state(i),-value)
      end do

      do i=1,self%nsurface_state
         _GET_HORIZONTAL_(self%id_surface_state(i),value)
         if (value/=i+surface_state_offset) call self%fatal_error('do_surface','invalid value of surface state variable.')
         _SET_SURFACE_ODE_(self%id_surface_state(i),-value)
      end do

      do i=1,self%nbottom_state
         _GET_HORIZONTAL_(self%id_bottom_state(i),value)
         if (value/=i+bottom_state_offset) call self%fatal_error('do_surface','invalid value of bottom state variable.')
      end do

      _GET_(self%id_depth,value)
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
       if (value/=0) call self%fatal_error('do_surface','depth should be 0.')
#else
       if (value/=2) call self%fatal_error('do_surface','depth should be 2.')
#endif

      _GET_(self%id_dep,value)
      if (value/=1+interior_dependency_offset) call self%fatal_error('do_surface','invalid value of interior dependency #1.')
      _GET_HORIZONTAL_(self%id_hz_dep,value)
      if (value/=1+horizontal_dependency_offset) call self%fatal_error('do_surface','invalid value of horizontal dependency #1.')

   _HORIZONTAL_LOOP_END_
end subroutine do_surface

subroutine do_bottom(self,_ARGUMENTS_DO_SURFACE_)
   class (type_test_model),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   integer  :: i
   real(rk) :: value

   _HORIZONTAL_LOOP_BEGIN_
      do i=1,self%nstate
         _GET_(self%id_state(i),value)
         if (value/=i+interior_state_offset) call self%fatal_error('do_bottom','invalid value of interior state variable.')
         _SET_BOTTOM_EXCHANGE_(self%id_state(i),-value)
      end do

      do i=1,self%nsurface_state
         _GET_HORIZONTAL_(self%id_surface_state(i),value)
         if (value/=i+surface_state_offset) call self%fatal_error('do_bottom','invalid value of surface state variable.')
      end do

      do i=1,self%nbottom_state
         _GET_HORIZONTAL_(self%id_bottom_state(i),value)
         if (value/=i+bottom_state_offset) call self%fatal_error('do_bottom','invalid value of bottom state variable.')
         _SET_BOTTOM_ODE_(self%id_bottom_state(i),-value)
      end do

      _GET_(self%id_depth,value)
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
       if (value/=1) call self%fatal_error('do_bottom','depth should be 1.')
#else
       if (value/=2) call self%fatal_error('do_bottom','depth should be 2.')
#endif

      _GET_(self%id_dep,value)
      if (value/=1+interior_dependency_offset) call self%fatal_error('do_bottom','invalid value of interior dependency #1.')
      _GET_HORIZONTAL_(self%id_hz_dep,value)
      if (value/=1+horizontal_dependency_offset) call self%fatal_error('do_bottom','invalid value of horizontal dependency #1.')

   _HORIZONTAL_LOOP_END_
end subroutine do_bottom

end module
   
program test_host

use fabm
use fabm_config
use fabm_driver

use test_models

implicit none

integer,parameter :: ntest = 1

#if _FABM_DIMENSION_COUNT_>0
integer :: _LOCATION_
#endif

#ifdef _INTERIOR_IS_VECTORIZED_
integer :: loop_start, loop_stop
#endif

#ifdef _HAS_MASK_
#  ifndef _FABM_HORIZONTAL_MASK_
_FABM_MASK_TYPE_,allocatable,target _DIMENSION_GLOBAL_ :: mask
real(rk),allocatable _DIMENSION_GLOBAL_ :: mask_values
#  endif
_FABM_MASK_TYPE_,allocatable,target _DIMENSION_GLOBAL_HORIZONTAL_ :: mask
real(rk),allocatable _DIMENSION_GLOBAL_HORIZONTAL_ :: mask_values
#endif

real(rk),allocatable,target _DIMENSION_GLOBAL_PLUS_1_            :: interior_state
real(rk),allocatable,target _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: surface_state
real(rk),allocatable,target _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: bottom_state

real(rk),allocatable _DIMENSION_EXT_SLICE_PLUS_1_        :: dy
real(rk),allocatable _DIMENSION_EXT_SLICE_PLUS_1_        :: w
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
i__=10000
#endif
#if _FABM_DIMENSION_COUNT_>1
j__=20
#endif
#if _FABM_DIMENSION_COUNT_>2
k__=12
#endif

#if _FABM_DIMENSION_COUNT_>0
domain_extent = (/ _LOCATION_ /)
#endif

#ifdef _INTERIOR_IS_VECTORIZED_
loop_start = 1
loop_stop = domain_extent(_FABM_VECTORIZED_DIMENSION_INDEX_)
#endif

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
#  ifdef _FABM_HORIZONTAL_MASK_
allocate(mask _INDEX_HORIZONTAL_LOCATION_)
allocate(mask_values _INDEX_HORIZONTAL_LOCATION_)
#  else
allocate(mask _INDEX_LOCATION_)
allocate(mask_values _INDEX_LOCATION_)
#  endif

call start_test('fabm_set_mask')
call fabm_set_mask(model,mask)
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
call model%set_bottom_index(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_))
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
      call random_number(mask_values)
      mask = _FABM_UNMASKED_VALUE_
      where (mask_values>0.5) mask = _FABM_UNMASKED_VALUE_+1
#endif
   end subroutine randomize_mask

   subroutine test_update
      call randomize_mask

      ! ======================================================================
      ! Initialize all state variables
      ! ======================================================================

      call start_test('fabm_initialize_state')
      do j__=1,domain_extent(2)
         call fabm_initialize_state(model _ARGUMENTS_INTERIOR_IN_)
      end do
      j__=domain_extent(2)
      call report_test_result()

      call start_test('fabm_initialize_bottom_state')
      call fabm_initialize_bottom_state(model _ARGUMENTS_HORIZONTAL_IN_)
      call report_test_result()

      call start_test('fabm_initialize_surface_state')
      call fabm_initialize_surface_state(model _ARGUMENTS_HORIZONTAL_IN_)
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
      do _VERTICAL_ITERATOR_=1,domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         depth(:,_VERTICAL_ITERATOR_) = real(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-_VERTICAL_ITERATOR_,rk)/(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-1)
#  else
         depth(:,_VERTICAL_ITERATOR_) = real(_VERTICAL_ITERATOR_-1,rk)/(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-1)
#  endif
      end do
      _VERTICAL_ITERATOR_=domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)
#else
      ! No depth dimension
      depth = 2
#endif
      call apply_mask_3d(depth,-999._rk-interior_dependency_offset)

      do ivar=1,size(model%state_variables)
         interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar) = ivar+interior_state_offset
         call apply_mask_3d(interior_state(:,:,ivar),-999._rk-interior_state_offset)
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
      dy = 0
      call fabm_do(model _ARGUMENTS_INTERIOR_IN_,dy)
      do ivar=1,size(model%state_variables)
         call check_interior_slice_plus_1(dy,ivar,0.0_rk,-real(ivar+interior_state_offset,rk))
      end do
      call report_test_result()

      ! ======================================================================
      ! Retrieve surface fluxes of interior state variables, source terms of surface-associated state variables.
      ! ======================================================================

      call start_test('fabm_do_surface')
      flux = 0
      sms_sf = 0
      call fabm_do_surface(model _ARGUMENTS_HORIZONTAL_IN_,flux,sms_sf)
      do ivar=1,size(model%state_variables)
         call check_horizontal_slice_plus_1(flux,ivar,0.0_rk,-real(ivar+interior_state_offset,rk))
      end do
      do ivar=1,size(model%surface_state_variables)
         call check_horizontal_slice_plus_1(sms_sf,ivar,0.0_rk,-real(ivar+surface_state_offset,rk))
      end do
      call report_test_result()

      ! ======================================================================
      ! Retrieve bottom fluxes of interior state variables, source terms of bottom-associated state variables.
      ! ======================================================================

      call start_test('fabm_do_bottom')
      flux = 0
      sms_bt = 0
      call fabm_do_bottom(model _ARGUMENTS_HORIZONTAL_IN_,flux,sms_bt)
      do ivar=1,size(model%state_variables)
         call check_horizontal_slice_plus_1(flux,ivar,0.0_rk,-real(ivar+interior_state_offset,rk))
      end do
      do ivar=1,size(model%bottom_state_variables)
         call check_horizontal_slice_plus_1(sms_bt,ivar,0.0_rk,-real(ivar+bottom_state_offset,rk))
      end do
      call report_test_result()

      ! ======================================================================
      ! Retrieve vertical velocities (sinking, floating, active movement).
      ! ======================================================================

      call start_test('fabm_get_vertical_movement')
      call fabm_get_vertical_movement(model _ARGUMENTS_INTERIOR_IN_,w)
      do ivar=1,size(model%state_variables)
         call check_interior_slice_plus_1(w,ivar,0.0_rk,-real(ivar+interior_state_offset,rk))
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

   subroutine apply_mask_3d(dat,missing_value)
      real(rk) _DIMENSION_GLOBAL_,intent(inout) :: dat
      real(rk),                   intent(in)    :: missing_value
#ifdef _FABM_HORIZONTAL_MASK_
      integer :: j__
      do j__=1,domain_extent(2)
         where (.not._IS_UNMASKED_(mask)) dat(:,j__) = missing_value
      end do
#else
      where (.not._IS_UNMASKED_(mask)) dat = missing_value
#endif
   end subroutine

   subroutine apply_mask_2d(dat,missing_value)
      real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,intent(inout) :: dat
      real(rk),                              intent(in)    :: missing_value
#ifdef _FABM_HORIZONTAL_MASK_
      where (.not._IS_UNMASKED_(mask)) dat = missing_value
#else
      where (.not._IS_UNMASKED_(mask)) dat = missing_value
#endif
   end subroutine

   subroutine check_interior_slice_plus_1(dat,index,required_masked_value,required_value)
      real(rk) _DIMENSION_EXT_SLICE_PLUS_1_,intent(in) :: dat
      integer,                              intent(in) :: index
      real(rk),                             intent(in) :: required_masked_value,required_value
#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_interior_slice(dat(:,ivar),required_masked_value,required_value)
#else
      call check_interior_slice(dat(ivar),required_masked_value,required_value)
#endif
   end subroutine

   subroutine check_interior_slice(dat,required_masked_value,required_value)
      real(rk) _DIMENSION_EXT_SLICE_,intent(in) :: dat
      real(rk),                      intent(in) :: required_masked_value,required_value
   end subroutine

   subroutine check_horizontal_slice_plus_1(dat,index,required_masked_value,required_value)
      real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(in) :: dat
      integer,                                     intent(in) :: index
      real(rk),                                    intent(in) :: required_masked_value,required_value
#ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_horizontal_slice(dat(:,ivar),required_masked_value,required_value)
#else
      call check_horizontal_slice(dat(ivar),required_masked_value,required_value)
#endif
   end subroutine

   subroutine check_horizontal_slice(dat,required_masked_value,required_value)
      real(rk) _DIMENSION_HORIZONTAL_SLICE_,intent(in) :: dat
      real(rk),                             intent(in) :: required_masked_value,required_value
#ifdef _HORIZONTAL_IS_VECTORIZED_
#  ifdef _HAS_MASK_
#    ifdef _FABM_HORIZONTAL_MASK_
      ! 1D slice; horizontal-only mask
      if (any(dat/=required_masked_value.and..not._IS_UNMASKED_(mask))) then
         call driver%fatal_error('check_horizontal_slice','one or more masked cells do not have the value required.')
      end if
      if (any(dat/=required_value.and._IS_UNMASKED_(mask))) then
         call driver%fatal_error('check_horizontal_slice','one or more non-masked cells do not have the value required.')
      end if
#    else
      ! 1D slice; interior mask [full domain]
      where (.not._IS_UNMASKED_(mask)) dat = missing_value
#    endif
#  else
      ! 1D slice; no mask
#  endif
#else
      ! 0D slice; no mask
      if (dat/=required_value) then
         call driver%fatal_error('check_horizontal_slice','variable does not have the value required.')
#endif
   end subroutine check_horizontal_slice

end program
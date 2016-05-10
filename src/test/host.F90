#include "fabm_driver.h"
#include "fabm_private.h"

module test_models

use fabm_types

type,extends(type_base_model) :: type_test_model
   type (type_state_variable_id) :: id_state
   type (type_surface_state_variable_id) :: id_surface_state
   type (type_dependency_id) :: id_dep
   type (type_horizontal_dependency_id) :: id_hz_dep
contains
   procedure :: initialize
   procedure :: do_surface
end type

   contains

subroutine initialize(self,configunit)
   class (type_test_model), intent(inout), target :: self
   integer,                 intent(in)            :: configunit
   call self%register_state_variable(self%id_state,'state','','state variable #1')
   call self%register_surface_state_variable(self%id_surface_state,'surface_state','','surface state variable #1')
   call self%register_dependency(self%id_dep,standard_variables%temperature)
   call self%register_dependency(self%id_hz_dep,standard_variables%wind_speed)
end subroutine initialize

subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_test_model),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: value

   _HORIZONTAL_LOOP_BEGIN_
      _GET_(self%id_state,value)
      if (value/=1) call self%fatal_error('do_surface','value of interior state variable #1 should be 1.')
      _GET_HORIZONTAL_(self%id_surface_state,value)
      if (value/=1001) call self%fatal_error('do_surface','value of surface state variable #1 should be 1001.')
      _GET_(self%id_dep,value)
      if (value/=3001) call self%fatal_error('do_surface','value of interior dependency #1 should be 3001.')
      _GET_HORIZONTAL_(self%id_hz_dep,value)
      if (value/=4001) call self%fatal_error('do_surface','value of horizontal dependency #1 should be 4001.')
   _HORIZONTAL_LOOP_END_
end subroutine do_surface

end module
   
program test_host

use fabm
use fabm_config

use test_models

implicit none

#if _FABM_DIMENSION_COUNT_>0
integer :: _LOCATION_
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

real(rk),allocatable _DIMENSION_EXT_SLICE_PLUS_1_ :: dy
real(rk),allocatable _DIMENSION_EXT_SLICE_PLUS_1_ :: w
real(rk),allocatable _DIMENSION_HORIZONTAL_SLICE_PLUS_1_ :: sms_sf,flux_sf
real(rk),allocatable _DIMENSION_HORIZONTAL_SLICE_PLUS_1_ :: sms_bt,flux_bt

real(rk),allocatable,target _DIMENSION_GLOBAL_ :: temp
real(rk),allocatable,target _DIMENSION_GLOBAL_HORIZONTAL_ :: wind_speed

type (type_model) :: model
class (type_test_model), pointer :: test_model
integer :: ivar

integer :: extents(_FABM_DIMENSION_COUNT_)

#if _FABM_DIMENSION_COUNT_>0
i__=10
#endif
#if _FABM_DIMENSION_COUNT_>1
j__=11
#endif
#if _FABM_DIMENSION_COUNT_>2
k__=12
#endif

extents = (/_LOCATION_/)

#if _FABM_DIMENSION_COUNT_>0
loop_start = 1
loop_stop = extents(_FABM_VECTORIZED_DIMENSION_INDEX_)
#endif

call fabm_initialize_library()

allocate(test_model)
call model%root%add_child(test_model,'test_model','test model',configunit=-1)

call fabm_initialize(model)

! Provide extents of the spatial domain.
call fabm_set_domain(model _ARGUMENTS_LOCATION_)

#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
allocate(mask _INDEX_HORIZONTAL_LOCATION_)
allocate(mask_values _INDEX_HORIZONTAL_LOCATION_)
#  else
allocate(mask _INDEX_LOCATION_)
allocate(mask_values _INDEX_LOCATION_)
#  endif
call random_number(mask_values)
mask = _FABM_UNMASKED_VALUE_
where (mask_values>0.5) mask = _FABM_UNMASKED_VALUE_+1
call fabm_set_mask(model,mask)
#endif

! Specify vertical index of surface and bottom
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
call model%set_bottom_index(1)
call model%set_surface_index(extents(_FABM_DEPTH_DIMENSION_INDEX_))
#  else
call model%set_surface_index(1)
call model%set_bottom_index(extents(_FABM_DEPTH_DIMENSION_INDEX_))
#  endif
#endif

allocate(interior_state(_PREARG_LOCATION_ size(model%state_variables)))
allocate(surface_state(_PREARG_HORIZONTAL_LOCATION_ size(model%surface_state_variables)))
allocate(bottom_state(_PREARG_HORIZONTAL_LOCATION_ size(model%bottom_state_variables)))

allocate(temp _INDEX_LOCATION_)
allocate(wind_speed _INDEX_HORIZONTAL_LOCATION_)

! Send pointers to state variable data to FABM
do ivar=1,size(model%state_variables)
   call model%link_interior_state_data(ivar,interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar))
end do
do ivar=1,size(model%surface_state_variables)
   call model%link_surface_state_data(ivar,surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar))
end do
do ivar=1,size(model%bottom_state_variables)
   call model%link_bottom_state_data(ivar,bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar))
end do

! Transfer pointer to environmental data
! Array temp with extents nx,ny,nz is assumed to be declared and accessible.
! Do this for all variables on FABM's standard variable list that the model can provide.
! For this list, visit http://fabm.net/standard_variables
call model%link_interior_data(standard_variables%temperature,temp)
call model%link_horizontal_data(standard_variables%wind_speed,wind_speed)

! Check whether FABM has all dependencies fulfilled
! (i.e., whether all required calls for fabm_link_*_data have been made)
call fabm_check_ready(model)

! Initialize the tracers
! This sets the values of arrays sent to fabm_link_bulk_state_data, in this case interior_state.
do j__=1,extents(2)
   call fabm_initialize_state(model _ARGUMENTS_INTERIOR_IN_)
end do
call fabm_initialize_bottom_state(model _ARGUMENTS_HORIZONTAL_IN_)
call fabm_initialize_surface_state(model _ARGUMENTS_HORIZONTAL_IN_)

j__=extents(2)

! At this point, initialization is complete.
! Routines below would typically be called every time step.
! This example only processes a single 1:nx slice.
! You would typically surround this logic with loops over all j and k.

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
allocate(dy(loop_start:loop_stop,size(model%state_variables)))
allocate(w(loop_start:loop_stop,size(model%state_variables)))
#else
allocate(dy(size(model%state_variables)))
allocate(w(size(model%state_variables)))
#endif

#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&(_FABM_DEPTH_DIMENSION_INDEX_!=_FABM_VECTORIZED_DIMENSION_INDEX_)
allocate(sms_sf(loop_start:loop_stop,size(model%surface_state_variables)))
allocate(flux_sf(loop_start:loop_stop,size(model%state_variables)))
allocate(sms_bt(loop_start:loop_stop,size(model%bottom_state_variables)))
allocate(flux_bt(loop_start:loop_stop,size(model%state_variables)))
#else
allocate(dy(size(model%state_variables)))
allocate(w(size(model%state_variables)))
#endif

temp = 3001
do j__=1,extents(2)
   where (.not._IS_UNMASKED_(mask)) temp(:,j__) = -3999
end do
j__=extents(2)

wind_speed = 4001
where (.not._IS_UNMASKED_(mask)) wind_speed = -4999

do ivar=1,size(model%state_variables)
   interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar) = ivar
   do j__=1,extents(2)
      where (.not._IS_UNMASKED_(mask)) interior_state(:,j__,ivar) = -999
   end do
   j__=extents(2)
end do
do ivar=1,size(model%surface_state_variables)
   surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = ivar+1000
   where (.not._IS_UNMASKED_(mask)) surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = -1999
end do
do ivar=1,size(model%bottom_state_variables)
   bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = ivar+2000
   where (.not._IS_UNMASKED_(mask)) bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = -2999
end do

! Retrieve tracer source terms (tracer units s-1).
! Array dy(1:nx,1:size(model%state_variables)) is assumed to be declared.
dy = 0
call fabm_do(model _ARGUMENTS_INTERIOR_IN_,dy)

flux_sf = 0
sms_sf = 0
call fabm_do_surface(model _ARGUMENTS_HORIZONTAL_IN_,flux_sf,sms_sf)

flux_bt = 0
sms_bt = 0
call fabm_do_bottom(model _ARGUMENTS_HORIZONTAL_IN_,flux_bt,sms_bt)

! Retrieve vertical velocities (sinking, floating, active movement) in m s-1.
! Array w(1:nx,1:size(model%state_variables)) is assumed to be declared.
call fabm_get_vertical_movement(model _ARGUMENTS_INTERIOR_IN_,w)

! Here you would time-integrate the advection-diffusion-reaction equations
! of all tracers, combining the transport terms with the biogeochemical source
! terms dy and vertical velocities w. This should result in an updated interior_state.
end program
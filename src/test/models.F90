#include "fabm_driver.h"

module test_models

use fabm_types

implicit none

private

integer, parameter, public :: interior_state_offset = 0
integer, parameter, public :: surface_state_offset  = 1000
integer, parameter, public :: bottom_state_offset   = 2000

integer, parameter, public :: interior_dependency_offset = 3000
integer, parameter, public :: horizontal_dependency_offset = 4000

real(rk), parameter :: epsilon = 1e-14_rk

type, extends(type_base_model), public :: type_test_model
   type (type_state_variable_id),        allocatable :: id_state(:)
   type (type_surface_state_variable_id),allocatable :: id_surface_state(:)
   type (type_bottom_state_variable_id), allocatable :: id_bottom_state(:)

   type (type_diagnostic_variable_id),           allocatable :: id_diag(:)
   type (type_horizontal_diagnostic_variable_id),allocatable :: id_horizontal_diag(:)

   type (type_dependency_id)            :: id_dep
   type (type_dependency_id)            :: id_depth
   type (type_horizontal_dependency_id) :: id_hz_dep

   integer :: nstate         = 12
   integer :: nsurface_state = 11
   integer :: nbottom_state  = 10
   integer :: nint_diag      = 12
   integer :: nsurface_diag  = 4
   integer :: nbottom_diag   = 6
   integer :: nhz_diag_vert  = 3
   integer :: nint_diag_vert = 2
contains
   procedure :: initialize
   procedure :: do
   procedure :: do_surface
   procedure :: do_bottom
   procedure :: do_column
   procedure :: get_vertical_movement
end type

integer, save, public :: interior_loop_count = 0
integer, save, public :: surface_loop_count = 0
integer, save, public :: bottom_loop_count = 0
integer, save, public :: column_loop_count = 0
integer, save, public :: vertical_movement_loop_count = 0

type (type_universal_standard_variable), parameter :: state_total = type_universal_standard_variable(name='state', conserved=.true., aggregate_variable=.true.)

   contains

subroutine initialize(self,configunit)
   class (type_test_model), intent(inout), target :: self
   integer,                 intent(in)            :: configunit

   integer          :: i
   character(len=8) :: strindex
   real(rk) :: w

   allocate(self%id_state(self%nstate))
   do i=1,self%nstate
      write (strindex,'(i0)') i
      w = 0
      if (mod(i, 2) /= 0) w = -real(i+interior_state_offset,rk)
      call self%register_state_variable(self%id_state(i),'state'//trim(strindex),'','state variable #'//trim(strindex),vertical_movement=w, initial_value=1._rk+i+interior_state_offset, missing_value=-999._rk-interior_state_offset-i)
      call self%add_to_aggregate_variable(state_total, self%id_state(i))
   end do
   allocate(self%id_surface_state(self%nsurface_state))
   do i=1,self%nsurface_state
      write (strindex,'(i0)') i
      call self%register_state_variable(self%id_surface_state(i),'surface_state'//trim(strindex),'','surface state variable #'//trim(strindex), initial_value=1._rk+i+surface_state_offset, missing_value=-999._rk-surface_state_offset-i)
      call self%add_to_aggregate_variable(state_total, self%id_surface_state(i))
   end do
   allocate(self%id_bottom_state(self%nbottom_state))
   do i=1,self%nbottom_state
      write (strindex,'(i0)') i
      call self%register_state_variable(self%id_bottom_state(i),'bottom_state'//trim(strindex),'','bottom state variable #'//trim(strindex), initial_value=1._rk+i+bottom_state_offset, missing_value=-999._rk-bottom_state_offset-i)
      call self%add_to_aggregate_variable(state_total, self%id_bottom_state(i))
   end do
   call self%register_dependency(self%id_dep,standard_variables%temperature)
   call self%register_dependency(self%id_depth,standard_variables%depth)
   call self%register_dependency(self%id_hz_dep,standard_variables%wind_speed)

   allocate(self%id_diag(self%nint_diag + self%nint_diag_vert))
   do i=1,self%nint_diag
      write (strindex,'(i0)') i
      call self%register_diagnostic_variable(self%id_diag(i),'diagnostic'//trim(strindex),'','diagnostic variable #'//trim(strindex),missing_value=-999._rk - i)
   end do
   do i=1,self%nint_diag_vert
      write (strindex,'(i0)') i
      call self%register_diagnostic_variable(self%id_diag(self%nint_diag + i),'vertical_diagnostic'//trim(strindex),'','vertical diagnostic variable #'//trim(strindex),missing_value=-1999._rk - i, source=source_do_column)
   end do
   allocate(self%id_horizontal_diag(self%nsurface_diag + self%nbottom_diag + self%nhz_diag_vert))
   do i=1,self%nsurface_diag
      write (strindex,'(i0)') i
      call self%register_diagnostic_variable(self%id_horizontal_diag(i),'surface_diagnostic'//trim(strindex),'','surface diagnostic variable #'//trim(strindex),missing_value=-2999._rk - i,source=source_do_surface)
   end do
   do i=1,self%nbottom_diag
      write (strindex,'(i0)') i
      call self%register_diagnostic_variable(self%id_horizontal_diag(self%nsurface_diag + i),'bottom_diagnostic'//trim(strindex),'','bottom diagnostic variable #'//trim(strindex),missing_value=-3999._rk - i,source=source_do_bottom)
   end do
   do i=1,self%nhz_diag_vert
      write (strindex,'(i0)') i
      call self%register_diagnostic_variable(self%id_horizontal_diag(self%nsurface_diag + self%nbottom_diag + i),'vert_hz_diagnostic'//trim(strindex),'','horizontal diagnostic variable set from do_column #'//trim(strindex),missing_value=-4999._rk - i,source=source_do_column)
   end do
end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)
   class (type_test_model),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

   integer  :: i
   real(rk) :: value

   _LOOP_BEGIN_
      do i=1,self%nstate
         _GET_(self%id_state(i),value)
         if (value/=i+interior_state_offset) call self%fatal_error('do','invalid value of interior state variable.')
         _SET_ODE_(self%id_state(i),-value)
      end do

      do i=1,self%nsurface_state
         _GET_HORIZONTAL_(self%id_surface_state(i),value)
         if (value/=i+surface_state_offset) call self%fatal_error('do','invalid value of surface state variable.')
      end do

      do i=1,self%nbottom_state
         _GET_HORIZONTAL_(self%id_bottom_state(i),value)
         if (value/=i+bottom_state_offset) call self%fatal_error('do','invalid value of bottom state variable.')
      end do

      _GET_(self%id_dep,value)
      if (value/=1+interior_dependency_offset) call self%fatal_error('do','invalid value of interior dependency #1.')
      _GET_HORIZONTAL_(self%id_hz_dep,value)
      if (value/=1+horizontal_dependency_offset) call self%fatal_error('do','invalid value of horizontal dependency #1.')

      do i=1,self%nint_diag
         _SET_DIAGNOSTIC_(self%id_diag(i),999._rk+i)
      end do

      interior_loop_count = interior_loop_count + 1
   _LOOP_END_
end subroutine do

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
       if (abs(value)>epsilon) call self%fatal_error('do_surface','depth should be 0.')
#else
       if (abs(value-2)>epsilon) call self%fatal_error('do_surface','depth should be 2.')
#endif

      _GET_(self%id_dep,value)
      if (value/=1+interior_dependency_offset) call self%fatal_error('do_surface','invalid value of interior dependency #1.')
      _GET_HORIZONTAL_(self%id_hz_dep,value)
      if (value/=1+horizontal_dependency_offset) call self%fatal_error('do_surface','invalid value of horizontal dependency #1.')

      do i=1,self%nsurface_diag
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_horizontal_diag(i),2999._rk+i)
      end do

      surface_loop_count = surface_loop_count + 1
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
       if (abs(value-1)>epsilon) call self%fatal_error('do_bottom','depth should be 1.')
#else
       if (abs(value-2)>epsilon) call self%fatal_error('do_bottom','depth should be 2.')
#endif

      _GET_(self%id_dep,value)
      if (value/=1+interior_dependency_offset) call self%fatal_error('do_bottom','invalid value of interior dependency #1.')
      _GET_HORIZONTAL_(self%id_hz_dep,value)
      if (value/=1+horizontal_dependency_offset) call self%fatal_error('do_bottom','invalid value of horizontal dependency #1.')

      do i=1,self%nbottom_diag
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_horizontal_diag(self%nsurface_diag + i),3999._rk+i)
      end do

      bottom_loop_count = bottom_loop_count + 1
   _HORIZONTAL_LOOP_END_
end subroutine do_bottom

subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
   class (type_test_model),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_COLUMN_

   integer  :: i
   real(rk) :: value
   real(rk) :: old_depth

   old_depth = -1._rk

   do i=1,self%nsurface_state
      _GET_HORIZONTAL_(self%id_surface_state(i),value)
      if (value/=i+surface_state_offset) call self%fatal_error('do_column','invalid value of surface state variable.')
   end do

   do i=1,self%nbottom_state
      _GET_HORIZONTAL_(self%id_bottom_state(i),value)
      if (value/=i+bottom_state_offset) call self%fatal_error('do_column','invalid value of bottom state variable.')
   end do

   _GET_HORIZONTAL_(self%id_hz_dep,value)
   if (value/=1+horizontal_dependency_offset) call self%fatal_error('do_column','invalid value of horizontal dependency #1.')

   _VERTICAL_LOOP_BEGIN_
      do i=1,self%nstate
         _GET_(self%id_state(i),value)
         if (value/=i+interior_state_offset) call self%fatal_error('do_column','invalid value of interior state variable.')
      end do

      _GET_(self%id_dep,value)
      if (value/=1+interior_dependency_offset) call self%fatal_error('do_column','invalid value of interior dependency #1.')

      _GET_(self%id_depth,value)
      if (value <= old_depth) &
          call self%fatal_error('do_column','depth is not increasing as expected.')
      old_depth = value

      do i=1,self%nint_diag_vert
         _SET_DIAGNOSTIC_(self%id_diag(self%nint_diag + i),1999._rk+i)
      end do

      column_loop_count = column_loop_count + 1
   _VERTICAL_LOOP_END_

   do i=1,self%nhz_diag_vert
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_horizontal_diag(self%nsurface_diag + self%nbottom_diag + i),4999._rk+i)
   end do
end subroutine do_column

subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
   class (type_test_model),intent(in) :: self
   _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

   integer  :: i
   real(rk) :: value

   _LOOP_BEGIN_
      do i=1,self%nstate
         _GET_(self%id_state(i),value)
         if (value/=i+interior_state_offset) call self%fatal_error('get_vertical_movement','invalid value of interior state variable.')
         if (mod(i, 2) == 0) _SET_VERTICAL_MOVEMENT_(self%id_state(i),real(i+interior_state_offset,rk))
      end do

      do i=1,self%nsurface_state
         _GET_HORIZONTAL_(self%id_surface_state(i),value)
         if (value/=i+surface_state_offset) call self%fatal_error('get_vertical_movement','invalid value of surface state variable.')
      end do

      do i=1,self%nbottom_state
         _GET_HORIZONTAL_(self%id_bottom_state(i),value)
         if (value/=i+bottom_state_offset) call self%fatal_error('get_vertical_movement','invalid value of bottom state variable.')
      end do

      _GET_(self%id_dep,value)
      if (value/=1+interior_dependency_offset) call self%fatal_error('get_vertical_movement','invalid value of interior dependency #1.')
      _GET_HORIZONTAL_(self%id_hz_dep,value)
      if (value/=1+horizontal_dependency_offset) call self%fatal_error('get_vertical_movement','invalid value of horizontal dependency #1.')

      vertical_movement_loop_count = vertical_movement_loop_count + 1
   _LOOP_END_
end subroutine get_vertical_movement

end module

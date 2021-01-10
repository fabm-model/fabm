#include "fabm_driver.h"
#include "fabm_private.h"

module fabm_c_integrate
#if _FABM_DIMENSION_COUNT_ == 0
   use iso_c_binding, only: c_ptr, c_f_pointer, c_loc, c_int, c_double

   use fabm_types, only: rke
   use fabm_c
   use fabm_driver, only: driver
   use fabm, only: status_start_done

   implicit none

contains

   subroutine integrate(pmodel, nt, ny, t_, y_ini_, y_, dt, do_surface, do_bottom, cell_thickness) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: integrate
      type (c_ptr),  value, intent(in) :: pmodel
      integer(c_int),value, intent(in) :: nt, ny
      real(c_double),target,intent(in) :: t_(*), y_ini_(*), y_(*)
      real(c_double),value, intent(in) :: dt
      integer(c_int),value, intent(in) :: do_surface, do_bottom
      real(c_double), target, intent(in) :: cell_thickness(*)

      type (type_model_wrapper), pointer :: model
      real(c_double), pointer :: t(:), y_ini(:), y(:,:)
      integer                :: it
      real(rke)              :: t_cur
      real(rke), target      :: y_cur(ny)
      real(rke)              :: dy(ny)
      logical                :: surface, bottom
      real(c_double) _ATTRIBUTES_GLOBAL_, pointer :: cell_thickness_

      call c_f_pointer(pmodel, model)
      if (model%p%status < status_start_done) then
         call driver%fatal_error('integrate', 'start has not been called yet.')
         return
      end if
      if (ny /= size(model%p%interior_state_variables) + size(model%p%surface_state_variables) &
         + size(model%p%bottom_state_variables)) then
         call driver%fatal_error('integrate', 'ny is wrong length')
         return
      end if

      call c_f_pointer(c_loc(t_), t, (/nt/))
      call c_f_pointer(c_loc(y_ini_), y_ini, (/ny/))
      call c_f_pointer(c_loc(y_), y, (/ny, nt/))

      surface = int2logical(do_surface)
      bottom = int2logical(do_bottom)
      cell_thickness_ => c_f_pointer_interior(model, cell_thickness)
      call model%p%link_all_interior_state_data(y_cur(1:size(model%p%interior_state_variables)))
      call model%p%link_all_surface_state_data(y_cur(size(model%p%interior_state_variables) + 1: &
         size(model%p%interior_state_variables) + size(model%p%surface_state_variables)))
      call model%p%link_all_bottom_state_data(y_cur(size(model%p%interior_state_variables) &
         + size(model%p%surface_state_variables) + 1:))

      it = 1
      t_cur = t(1)
      y_cur = y_ini
      do while (it <= nt)
          if (t_cur >= t(it)) then
              y(:, it) = y_cur
              it = it + 1
          end if

          call model%p%prepare_inputs(t_cur)
          dy = 0.0_rke
          if (surface) call model%p%get_surface_sources(dy(1:size(model%p%interior_state_variables)), &
             dy(size(model%p%interior_state_variables) + 1:size(model%p%interior_state_variables) &
             + size(model%p%surface_state_variables)))
          if (bottom) call model%p%get_bottom_sources(dy(1:size(model%p%interior_state_variables)), &
             dy(size(model%p%interior_state_variables) + size(model%p%surface_state_variables) + 1:))
          if (surface .or. bottom) dy(1:size(model%p%interior_state_variables)) = dy(1:size(model%p%interior_state_variables)) &
             / cell_thickness_
          call model%p%get_interior_sources(dy(1:size(model%p%interior_state_variables)))
          y_cur = y_cur + dt * dy * 86400
          t_cur = t_cur + dt
      end do
   end subroutine integrate
#endif
end module

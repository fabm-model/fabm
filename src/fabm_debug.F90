#include "fabm_driver.h"
#include "fabm_private.h"

module fabm_debug

   use fabm_driver
   use fabm_types, only: rke

   implicit none

   private

   public check_interior_location, check_horizontal_location, check_vertical_location
   public check_extents_1d, check_extents_2d, check_extents_3d

contains

   subroutine check_interior_location(domain_start, domain_stop _POSTARG_INTERIOR_IN_, routine)
      integer, dimension(_FABM_DIMENSION_COUNT_), intent(in) :: domain_start, domain_stop
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      character(len=*), intent(in) :: routine

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_loop(_START_, _STOP_, domain_start(_FABM_VECTORIZED_DIMENSION_INDEX_), domain_stop(_FABM_VECTORIZED_DIMENSION_INDEX_), routine)
#endif
#if _FABM_DIMENSION_COUNT_>0&&_FABM_VECTORIZED_DIMENSION_INDEX_!=1
      call check_index(i__, domain_start(1), domain_stop(1), routine, 'i')
#endif
#if _FABM_DIMENSION_COUNT_>1&&_FABM_VECTORIZED_DIMENSION_INDEX_!=2
      call check_index(j__, domain_start(2), domain_stop(2), routine, 'j')
#endif
#if _FABM_DIMENSION_COUNT_>2&&_FABM_VECTORIZED_DIMENSION_INDEX_!=3
      call check_index(k__, domain_start(3), domain_stop(3), routine, 'k')
#endif
   end subroutine check_interior_location

   subroutine check_horizontal_location(domain_start, domain_stop _POSTARG_HORIZONTAL_IN_, routine)
      integer, dimension(_FABM_DIMENSION_COUNT_), intent(in) :: domain_start, domain_stop
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      character(len=*), intent(in) :: routine

#ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_loop(_START_, _STOP_, domain_start(_FABM_VECTORIZED_DIMENSION_INDEX_), domain_stop(_FABM_VECTORIZED_DIMENSION_INDEX_), routine)
#endif
#if _FABM_DIMENSION_COUNT_>0&&_FABM_VECTORIZED_DIMENSION_INDEX_!=1&&_FABM_DEPTH_DIMENSION_INDEX_!=1
      call check_index(i__, domain_start(1), domain_stop(1), routine, 'i')
#endif
#if _FABM_DIMENSION_COUNT_>1&&_FABM_VECTORIZED_DIMENSION_INDEX_!=2&&_FABM_DEPTH_DIMENSION_INDEX_!=2
      call check_index(j__, domain_start(2), domain_stop(2), routine, 'j')
#endif
#if _FABM_DIMENSION_COUNT_>2&&_FABM_VECTORIZED_DIMENSION_INDEX_!=3&&_FABM_DEPTH_DIMENSION_INDEX_!=3
      call check_index(k__, domain_start(3), domain_stop(3), routine, 'k')
#endif
   end subroutine check_horizontal_location

   subroutine check_vertical_location(domain_start, domain_stop _POSTARG_VERTICAL_IN_, routine)
      integer, dimension(_FABM_DIMENSION_COUNT_), intent(in) :: domain_start, domain_stop
      _DECLARE_ARGUMENTS_VERTICAL_IN_
      character(len=*), intent(in) :: routine

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      call check_loop(_VERTICAL_START_, _VERTICAL_STOP_, domain_start(_FABM_DEPTH_DIMENSION_INDEX_), domain_stop(_FABM_DEPTH_DIMENSION_INDEX_), routine)
#endif
#if _FABM_DIMENSION_COUNT_>0&&_FABM_DEPTH_DIMENSION_INDEX_!=1
      call check_index(i__, domain_start(1), domain_stop(1), routine, 'i')
#endif
#if _FABM_DIMENSION_COUNT_>1&&_FABM_DEPTH_DIMENSION_INDEX_!=2
      call check_index(j__, domain_start(2), domain_stop(2), routine, 'j')
#endif
#if _FABM_DIMENSION_COUNT_>2&&_FABM_DEPTH_DIMENSION_INDEX_!=3
      call check_index(k__, domain_start(3), domain_stop(3), routine, 'k')
#endif
   end subroutine check_vertical_location

   subroutine check_extents_1d(array, required_size1, routine, array_name, shape_description)
      real(rke),        intent(in) :: array(:)
      integer,          intent(in) :: required_size1
      character(len=*), intent(in) :: routine, array_name, shape_description

      character(len=8) :: actual, required

      if (size(array,1) /= required_size1) then
         write (actual,  '(i0)') size(array, 1)
         write (required,'(i0)') required_size1
         call fatal_error(routine, 'shape of argument ' // trim(array_name) // ' is (' // trim(actual) // ') but should be (' // trim(required) // ') = ' // trim(shape_description))
      end if
   end subroutine check_extents_1d

   subroutine check_extents_2d(array, required_size1, required_size2, routine, array_name, shape_description)
      real(rke),        intent(in) :: array(:,:)
      integer,          intent(in) :: required_size1, required_size2
      character(len=*), intent(in) :: routine, array_name, shape_description

      character(len=17) :: actual, required

      if (size(array,1) /= required_size1 .or. size(array,2) /= required_size2) then
         write (actual,  '(i0,a,i0)') size(array,1), ',', size(array,2)
         write (required,'(i0,a,i0)') required_size1, ',', required_size2
         call fatal_error(routine, 'shape of argument ' // trim(array_name) // ' is (' // trim(actual) // ') but should be (' // trim(required) // ') = ' // trim(shape_description))
      end if
   end subroutine check_extents_2d

   subroutine check_extents_3d(array, required_size1, required_size2, required_size3, routine, array_name, shape_description)
      real(rke),       intent(in) :: array(:,:,:)
      integer,         intent(in) :: required_size1, required_size2, required_size3
      character(len=*),intent(in) :: routine, array_name, shape_description

      character(len=26) :: actual, required

      if (size(array,1) /= required_size1 .or. size(array,2) /= required_size2 .or. size(array,3) /= required_size3) then
         write (actual,  '(i0,a,i0,a,i0)') size(array,1), ',', size(array,2), ',', size(array,3)
         write (required,'(i0,a,i0,a,i0)') required_size1, ',', required_size2, ',', required_size3
         call fatal_error(routine, 'shape of argument ' // trim(array_name) // ' is (' // trim(actual) // ') but should be (' // trim(required) // ') = ' // trim(shape_description))
      end if
   end subroutine check_extents_3d

   subroutine check_loop(istart, istop, imin, imax, routine)
      integer,          intent(in) :: istart, istop, imin, imax
      character(len=*), intent(in) :: routine

      character(len=8) :: str1, str2

      if (istart < imin) then
         write (str1,'(i0)') istart
         write (str2,'(i0)') imin
         call fatal_error(routine, 'Loop start index ' // trim(str1) // ' lies below minimum of ' // trim(str2) // '.')
      end if
      if (istop > imax) then
         write (str1,'(i0)') istop
         write (str2,'(i0)') imax
         call fatal_error(routine, 'Loop stop index ' // trim(str1) // ' exceeds maximum of ' // trim(str2) // '.')
      end if
      if (istart > istop) then
         write (str1,'(i0)') istart
         write (str2,'(i0)') istop
         call fatal_error(routine, 'Loop start index ' // trim(str1) // ' exceeds stop index ' // trim(str2) // '.')
      end if
   end subroutine check_loop

   subroutine check_index(i, imin, imax, routine, name)
      integer,          intent(in) :: i, imin, imax
      character(len=*), intent(in) :: routine, name

      character(len=8) :: str1, str2

      if (i < imin) then
         write (str1,'(i0)') i
         write (str2,'(i0)') imin
         call fatal_error(routine, 'Index ' // name // ' = ' // trim(str1) // ' lies below minimum of ' // trim(str2) // '.')
      end if
      if (i > imax) then
         write (str1,'(i0)') i
         write (str2,'(i0)') imax
         call fatal_error(routine, 'Index ' // name // ' = ' // trim(str1) // ' exceeds maximum of ' // trim(str2) // '.')
      end if
   end subroutine check_index

end module

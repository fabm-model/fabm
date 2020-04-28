#include "fabm_driver.h"
#include "fabm_private.h"

! This module define standard expressions that can be used by biogeochemical models.

module fabm_expressions

   use fabm_types
   use fabm_driver

   implicit none

   private

   public temporal_mean, vertical_mean, vertical_integral
   public type_interior_temporal_mean, type_horizontal_temporal_mean, type_vertical_integral

   type, extends(type_interior_expression) :: type_interior_temporal_mean
      real(rk) :: period   ! Time period to average over (s)
      integer  :: n
      real(rk) :: last_time, next_save_time
      integer  :: ioldest = -1
      real(rk) :: missing_value = -2.e20_rk
      logical  :: valid = .false.

      type (type_link), pointer :: link => null()
      integer :: in = -1
      real(rke),allocatable _DIMENSION_GLOBAL_PLUS_1_ :: history
   end type

   type, extends(type_horizontal_expression) :: type_horizontal_temporal_mean
      real(rk) :: period   ! Time period to average over (s)
      integer  :: n
      real(rk) :: last_time, next_save_time
      integer  :: ioldest = -1

      type (type_link), pointer :: link => null()
      integer :: in = -1
      real(rke),allocatable _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: history
   end type

   type, extends(type_horizontal_expression) :: type_vertical_integral
      real(rk) :: minimum_depth = 0.0_rk        ! Depth below surface in m (positive)
      real(rk) :: maximum_depth = huge(1.0_rk)  ! Depth below surface in m (positive)
      logical  :: average       = .false.       ! Whether to divide the depth integral by water depth, thus computing the vertical average
      character(len=attribute_length) :: input_name = ''

      type (type_link), pointer :: link => null()
   end type

   interface temporal_mean
      module procedure interior_temporal_mean
      module procedure horizontal_temporal_mean
   end interface

   interface vertical_mean
      module procedure vertical_dependency_mean
      module procedure vertical_state_mean
   end interface

   interface vertical_integral
      module procedure vertical_dependency_integral
      module procedure vertical_state_integral
   end interface

contains

   function vertical_dependency_mean(input, minimum_depth, maximum_depth) result(expression)
      type (type_dependency_id), intent(inout),target   :: input
      real(rk),                  intent(in),   optional :: minimum_depth,maximum_depth
      type (type_vertical_integral)                     :: expression
      expression = vertical_integral(input,minimum_depth,maximum_depth,average=.true.)
   end function

   function vertical_state_mean(input, minimum_depth, maximum_depth) result(expression)
      type (type_state_variable_id), intent(inout),target   :: input
      real(rk),                      intent(in),   optional :: minimum_depth, maximum_depth
      type (type_vertical_integral)                         :: expression
      expression = vertical_integral_generic(input,minimum_depth,maximum_depth,average=.true.)
   end function

   function vertical_dependency_integral(input, minimum_depth, maximum_depth, average) result(expression)
      type (type_dependency_id), intent(inout),target   :: input
      real(rk),                  intent(in),   optional :: minimum_depth, maximum_depth
      logical,                   intent(in),   optional :: average
      type (type_vertical_integral)                     :: expression
      expression = vertical_integral_generic(input, minimum_depth, maximum_depth, average)
   end function

   function vertical_state_integral(input, minimum_depth, maximum_depth, average) result(expression)
      type (type_state_variable_id), intent(inout),target   :: input
      real(rk),                      intent(in),   optional :: minimum_depth, maximum_depth
      logical,                       intent(in),   optional :: average
      type (type_vertical_integral)                         :: expression
      expression = vertical_integral_generic(input, minimum_depth, maximum_depth, average)
   end function

   function vertical_integral_generic(input, minimum_depth, maximum_depth, average) result(expression)
      class (type_variable_id), intent(inout),target   :: input
      real(rk),                 intent(in),   optional :: minimum_depth, maximum_depth
      logical,                  intent(in),   optional :: average
      type (type_vertical_integral)                    :: expression

      character(len=attribute_length) :: postfix

      if (.not.associated(input%link)) call fatal_error('fabm_expressions::vertical_mean', &
         'Input variable has not been registered yet.')
      expression%input_name = input%link%target%name

      ! Create a name for the expression
      postfix = ''
      if (present(minimum_depth).and.present(maximum_depth)) then
         if (minimum_depth>maximum_depth) call fatal_error('fabm_expressions::vertical_mean', &
            'Minimum depth exceeds maximum depth.')
         write (postfix,'(a,i0,a,i0,a)') '_between_',int(minimum_depth),'_m_and_',int(maximum_depth),'_m'
      elseif (present(minimum_depth)) then
         write (postfix,'(a,i0,a)') '_below_',int(minimum_depth),'_m'
      elseif (present(maximum_depth)) then
         write (postfix,'(a,i0,a)') '_above_',int(maximum_depth),'_m'
      end if
      if (present(average)) expression%average = average

      if (expression%average) then
         expression%output_name = 'vertical_mean_'//trim(input%link%name)//trim(postfix)
      else
         expression%output_name = 'integral_of_'//trim(input%link%name)//'_wrt_depth'//trim(postfix)
      end if

      expression%link => input%link
      if (present(minimum_depth)) expression%minimum_depth = minimum_depth
      if (present(maximum_depth)) expression%maximum_depth = maximum_depth
   end function

   function interior_temporal_mean(input, period, resolution, missing_value) result(expression)
      type (type_dependency_id), intent(inout), target :: input
      real(rk),                  intent(in)            :: period, resolution
      real(rk),optional,         intent(in)            :: missing_value
      type (type_interior_temporal_mean)               :: expression

      character(len=attribute_length) :: prefix, postfix

      if (.not.associated(input%link)) call fatal_error('fabm_expressions::interior_temporal_mean', &
         'Input variable has not been registered yet.')

      ! Create a name for the expression
      write (prefix,'(i0,a)') int(period),'_s_mean_'
      write (postfix,'(a,i0,a)') '_at_',int(resolution),'_s_resolution'
      expression%output_name = trim(prefix)//trim(input%link%name)//trim(postfix)

      expression%link => input%link
      expression%n = nint(period/resolution)
      expression%period = period
      if (present(missing_value)) expression%missing_value = missing_value
   end function

   function horizontal_temporal_mean(input, period, resolution) result(expression)
      type (type_horizontal_dependency_id), intent(inout), target :: input
      real(rk),                             intent(in)            :: period, resolution
      type (type_horizontal_temporal_mean)                        :: expression

      character(len=attribute_length) :: prefix, postfix

      if (.not.associated(input%link)) call fatal_error('fabm_expressions::horizontal_temporal_mean', &
         'Input variable has not been registered yet.')

      ! Create a name for the expression
      write (prefix,'(i0,a)') int(period),'_s_mean_'
      write (postfix,'(a,i0,a)') '_at_',int(resolution),'_s_resolution'
      expression%output_name = trim(prefix)//trim(input%link%name)//trim(postfix)

      expression%link => input%link
      expression%n = nint(period/resolution)
      expression%period = period
   end function

end module fabm_expressions

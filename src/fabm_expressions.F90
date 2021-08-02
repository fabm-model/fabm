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
      real(rk) :: missing_value = -2.e20_rk
      type (type_link), pointer :: link => null()
      integer :: in = -1

      real(rk), private :: previous_time, bin_end_time
      integer,  private :: ioldest  = -1
      integer,  private :: icurrent = -1
      logical,  private :: complete = .false.
      real(rke), allocatable _DIMENSION_GLOBAL_PLUS_1_ :: history
#if _FABM_DIMENSION_COUNT_>0
      real(rke), allocatable _DIMENSION_GLOBAL_ :: previous_value, last_exact_mean, mean
#else   
      real(rke) :: previous_value, last_exact_mean, mean
#endif
   contains
      procedure :: update => interior_temporal_mean_update
   end type

   type, extends(type_horizontal_expression) :: type_horizontal_temporal_mean
      real(rk) :: period   ! Time period to average over (s)
      integer  :: n
      real(rk) :: last_time, next_save_time
      integer  :: ioldest = -1
      integer  :: icurrent = -1

      type (type_link), pointer :: link => null()
      integer :: in = -1
      real(rke), allocatable _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: history
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
      type (type_dependency_id), intent(inout), target :: input
      real(rk), optional,        intent(in)            :: minimum_depth,maximum_depth
      type (type_vertical_integral)                    :: expression
      expression = vertical_integral(input,minimum_depth,maximum_depth,average=.true.)
   end function

   function vertical_state_mean(input, minimum_depth, maximum_depth) result(expression)
      type (type_state_variable_id), intent(inout), target :: input
      real(rk), optional,            intent(in)            :: minimum_depth, maximum_depth
      type (type_vertical_integral)                        :: expression
      expression = vertical_integral_generic(input,minimum_depth,maximum_depth,average=.true.)
   end function

   function vertical_dependency_integral(input, minimum_depth, maximum_depth, average) result(expression)
      type (type_dependency_id), intent(inout),target :: input
      real(rk), optional,        intent(in)           :: minimum_depth, maximum_depth
      logical,  optional,        intent(in)           :: average
      type (type_vertical_integral)                     :: expression
      expression = vertical_integral_generic(input, minimum_depth, maximum_depth, average)
   end function

   function vertical_state_integral(input, minimum_depth, maximum_depth, average) result(expression)
      type (type_state_variable_id), intent(inout), target :: input
      real(rk), optional,            intent(in)            :: minimum_depth, maximum_depth
      logical,  optional,            intent(in)            :: average
      type (type_vertical_integral)                          :: expression
      expression = vertical_integral_generic(input, minimum_depth, maximum_depth, average)
   end function

   function vertical_integral_generic(input, minimum_depth, maximum_depth, average) result(expression)
      class (type_variable_id), intent(inout), target :: input
      real(rk), optional,       intent(in)            :: minimum_depth, maximum_depth
      logical,  optional,       intent(in)            :: average
      type (type_vertical_integral)                   :: expression

      character(len=attribute_length) :: postfix

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::vertical_mean', &
         'Input variable has not been registered yet.')
      expression%input_name = input%link%target%name

      ! Create a name for the expression
      postfix = ''
      if (present(minimum_depth) .and. present(maximum_depth)) then
         if (minimum_depth>maximum_depth) call fatal_error('fabm_expressions::vertical_mean', &
            'Minimum depth exceeds maximum depth.')
         write (postfix,'(a,i0,a,i0,a)') '_between_', int(minimum_depth), '_m_and_', int(maximum_depth), '_m'
      elseif (present(minimum_depth)) then
         write (postfix,'(a,i0,a)') '_below_', int(minimum_depth), '_m'
      elseif (present(maximum_depth)) then
         write (postfix,'(a,i0,a)') '_above_', int(maximum_depth), '_m'
      end if
      if (present(average)) expression%average = average

      if (expression%average) then
         expression%output_name = 'vertical_mean_' // trim(input%link%name) // trim(postfix)
      else
         expression%output_name = 'integral_of_' // trim(input%link%name) // '_wrt_depth' // trim(postfix)
      end if

      expression%link => input%link
      if (present(minimum_depth)) expression%minimum_depth = minimum_depth
      if (present(maximum_depth)) expression%maximum_depth = maximum_depth
   end function

   function interior_temporal_mean(input, period, resolution, missing_value) result(expression)
      type (type_dependency_id), intent(inout), target :: input
      real(rk),                  intent(in)            :: period, resolution
      real(rk), optional,        intent(in)            :: missing_value
      type (type_interior_temporal_mean)               :: expression

      character(len=attribute_length) :: prefix, postfix

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::interior_temporal_mean', &
         'Input variable has not been registered yet.')

      ! Create a name for the expression
      write (prefix,'(i0,a)') int(period), '_s_mean_'
      write (postfix,'(a,i0,a)') '_at_', int(resolution), '_s_resolution'
      expression%output_name = trim(prefix) // trim(input%link%name) // trim(postfix)

      expression%link => input%link
      expression%n = nint(period / resolution)
      expression%period = period
      if (present(missing_value)) expression%missing_value = missing_value
   end function

   function horizontal_temporal_mean(input, period, resolution) result(expression)
      class (type_horizontal_dependency_id), intent(inout), target :: input
      real(rk),                              intent(in)            :: period, resolution
      type (type_horizontal_temporal_mean)                         :: expression

      character(len=attribute_length) :: prefix, postfix

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::horizontal_temporal_mean', &
         'Input variable has not been registered yet.')

      ! Create a name for the expression
      write (prefix,'(i0,a)') int(period), '_s_mean_'
      write (postfix,'(a,i0,a)') '_at_', int(resolution), '_s_resolution'
      expression%output_name = trim(prefix) // trim(input%link%name) // trim(postfix)

      expression%link => input%link
      expression%n = nint(period / resolution)
      expression%period = period
   end function

   subroutine interior_temporal_mean_update(self, time, value _POSTARG_LOCATION_RANGE_)
      class (type_interior_temporal_mean), intent(inout) :: self
      real(rke),                           intent(in)    :: time
      real(rke) _ATTRIBUTES_GLOBAL_,       intent(in)    :: value
      _DECLARE_ARGUMENTS_LOCATION_RANGE_

      integer  :: ibin
      real(rke) :: dt, w, dt_bin
      _DECLARE_LOCATION_

      ! Note that all array processing below uses explicit loops in order to respect
      ! any limits on the active domain given by the _LOCATION_RANGE_ argument.

      dt_bin = self%period / self%n

      if (self%ioldest == -1) then
         ! Start of simulation
         self%previous_time = time
         self%bin_end_time = time + dt_bin
         self%icurrent = 1
         self%ioldest = 2
      end if

      do while (time >= self%bin_end_time)
         ! Linearly interpolate to value at end-of-bin time and add that to the current bin
         dt = self%bin_end_time - self%previous_time
         w = dt / (time - self%previous_time)   ! weight for current time (leaving 1-w for previous time)
         _BEGIN_GLOBAL_LOOP_
            self%history(_PREARG_LOCATION_ self%icurrent) = self%history(_PREARG_LOCATION_ self%icurrent) &
               + ((1._rke - 0.5_rke * w) * self%previous_value _INDEX_LOCATION_ + 0.5_rke * w * value _INDEX_LOCATION_) &
               * dt / self%period
         _END_GLOBAL_LOOP_

         if (self%complete) then
            ! We already had a complete history (bins covering the full window size). Add the newly full bin, subtract the oldest bin
            _BEGIN_GLOBAL_LOOP_
               self%last_exact_mean _INDEX_LOCATION_ = self%last_exact_mean _INDEX_LOCATION_ &
                  - self%history(_PREARG_LOCATION_ self%ioldest) + self%history(_PREARG_LOCATION_ self%icurrent)
            _END_GLOBAL_LOOP_
         elseif (self%icurrent == self%n) then
            ! We just completed our history. Create the mean by summing all filled bins.
            do ibin = 1, self%n
               _BEGIN_GLOBAL_LOOP_
                  self%last_exact_mean _INDEX_LOCATION_ = self%last_exact_mean _INDEX_LOCATION_  &
                     + self%history(_PREARG_LOCATION_ ibin)
               _END_GLOBAL_LOOP_
            end do
            self%complete = .true.
         end if

         ! Update previous time and value to match end of current bin
         self%previous_time = self%bin_end_time
         _BEGIN_GLOBAL_LOOP_
            self%previous_value _INDEX_LOCATION_ = (1._rke - w) * self%previous_value _INDEX_LOCATION_ + w * value _INDEX_LOCATION_
         _END_GLOBAL_LOOP_

         ! Move to next bin: update indices, end time and empty newly current bin
         self%icurrent = self%ioldest
         self%ioldest = self%ioldest + 1
         if (self%ioldest > self%n + 1) self%ioldest = 1
         self%bin_end_time = self%bin_end_time + dt_bin
         _BEGIN_GLOBAL_LOOP_
            self%history(_PREARG_LOCATION_ self%icurrent) = 0
         _END_GLOBAL_LOOP_
      end do

      ! Compute average of previous and current value, multiply by time difference, pre-divide by window size, and add to current bin.
      _BEGIN_GLOBAL_LOOP_
         self%history(_PREARG_LOCATION_ self%icurrent) = self%history(_PREARG_LOCATION_ self%icurrent) &
            + 0.5_rke * (self%previous_value _INDEX_LOCATION_ + value _INDEX_LOCATION_) &
            * (time - self%previous_time) / self%period
      _END_GLOBAL_LOOP_

      if (self%complete) then
         ! We have a full history covering at least one window size. Update the running mean.
         ! The result is an approximation that assumes linear change over the period covered by the oldest bin.
         _BEGIN_GLOBAL_LOOP_
            self%mean _INDEX_LOCATION_ = self%last_exact_mean _INDEX_LOCATION_ &
               + self%history(_PREARG_LOCATION_ self%icurrent) &
               - self%history(_PREARG_LOCATION_ self%ioldest) * (time - self%bin_end_time + dt_bin) / dt_bin
         _END_GLOBAL_LOOP_
      end if

      ! Store current value to enable linear interpolation in subsequent call.
      self%previous_time = time
      _BEGIN_GLOBAL_LOOP_
         self%previous_value _INDEX_LOCATION_ = value _INDEX_LOCATION_
      _END_GLOBAL_LOOP_
   end subroutine

end module fabm_expressions

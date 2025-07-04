#include "fabm_driver.h"
#include "fabm_private.h"

! This module define standard expressions that can be used by biogeochemical models.

module fabm_expressions

   use fabm_types
   use fabm_driver
   use fabm_builtin_depth_integral

   implicit none

   private

   public temporal_mean, temporal_maximum, vertical_mean, vertical_integral
   public type_interior_temporal_mean, type_horizontal_temporal_mean, type_horizontal_temporal_maximum, type_vertical_integral

   type, extends(type_variable_id) :: type_global_interior_variable_id
      real(rke), pointer _ATTRIBUTES_GLOBAL_CONTIGUOUS_ :: p => null()
   end type

   type, extends(type_variable_id) :: type_global_horizontal_variable_id
      real(rke), pointer _ATTRIBUTES_GLOBAL_HORIZONTAL_CONTIGUOUS_ :: p => null()
   end type

   type, extends(type_variable_id) :: type_global_scalar_variable_id
      real(rke), pointer :: p => null()
   end type

   type, extends(type_interior_expression) :: type_interior_temporal_mean
      real(rk) :: period   ! Time period to average over (s)
      integer  :: n
      real(rk) :: missing_value = -2.e20_rk
      logical  :: use_incomplete_result = .false.

      type (type_link), pointer :: link => null()
      integer :: in = -1

      type (type_global_interior_variable_id), private :: previous_value, last_exact_mean, mean
      type (type_global_interior_variable_id), allocatable, private :: history(:)
      type (type_global_scalar_variable_id), private :: previous_time, start_time, icurrent
   contains
      procedure :: initialize => interior_temporal_mean_initialize
      procedure :: set_data   => interior_temporal_mean_set_data
      procedure :: update     => interior_temporal_mean_update
   end type

   type, extends(type_horizontal_expression) :: type_horizontal_temporal_mean
      real(rk) :: period   ! Time period to average over (s)
      integer  :: n
      real(rk) :: missing_value = 0.0_rk
      logical  :: use_incomplete_result = .false.

      type (type_link), pointer :: link => null()
      integer :: in = -1

      type (type_global_horizontal_variable_id), private :: previous_value, last_exact_mean, mean
      type (type_global_horizontal_variable_id), allocatable, private :: history(:)
      type (type_global_scalar_variable_id), private :: previous_time, start_time, icurrent
   contains
      procedure :: initialize => horizontal_temporal_mean_initialize
      procedure :: set_data   => horizontal_temporal_mean_set_data
      procedure :: update     => horizontal_temporal_mean_update
   end type

   type, extends(type_horizontal_expression) :: type_horizontal_temporal_maximum
      real(rk) :: period                     ! Time window to compute running maximum over (s)
      integer  :: n                          ! Number of bins to use to cover the period
      real(rk) :: missing_value = -2.e20_rk  ! Missing value to use until the simulation has covered the window size [period]

      type (type_link), pointer :: link => null()
      integer :: in = -1

      type (type_global_horizontal_variable_id), private :: previous_value, maximum
      type (type_global_horizontal_variable_id), allocatable, private :: history(:)
      type (type_global_scalar_variable_id), private :: previous_time, start_time, current
   contains
      procedure :: initialize => horizontal_temporal_maximum_initialize
      procedure :: set_data   => horizontal_temporal_maximum_set_data
      procedure :: update     => horizontal_temporal_maximum_update
   end type

   type, extends(type_horizontal_expression) :: type_vertical_integral
      real(rk) :: minimum_depth = 0.0_rk        ! Depth below surface in m (positive)
      real(rk) :: maximum_depth = huge(1.0_rk)  ! Depth below surface in m (positive)
      logical  :: average       = .false.       ! Whether to divide the depth integral by water depth, thus computing the vertical average

      type (type_link), pointer :: link => null()
   contains
      procedure :: initialize => vertical_integral_initialize
   end type

   interface temporal_mean
      module procedure interior_temporal_mean
      module procedure horizontal_temporal_mean
   end interface

   interface temporal_maximum
      module procedure horizontal_temporal_maximum
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
      expression = vertical_integral(input, minimum_depth, maximum_depth, average=.true.)
   end function

   function vertical_state_mean(input, minimum_depth, maximum_depth) result(expression)
      type (type_state_variable_id), intent(inout), target :: input
      real(rk), optional,            intent(in)            :: minimum_depth, maximum_depth
      type (type_vertical_integral)                        :: expression
      expression = vertical_integral_generic(input, minimum_depth, maximum_depth, average=.true.)
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

      ! Create a name for the expression
      postfix = ''
      if (present(minimum_depth) .and. present(maximum_depth)) then
         if (minimum_depth > maximum_depth) call fatal_error('fabm_expressions::vertical_mean', &
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

   function vertical_integral_initialize(self, model) result(link)
      class (type_vertical_integral), intent(inout) :: self
      class (type_base_model),        intent(inout) :: model
      type (type_link), pointer :: link

      class (type_depth_integral),         pointer :: integral
      class (type_bounded_depth_integral), pointer :: bounded_integral

      if (self%minimum_depth <= 0._rk .and. self%maximum_depth == huge(self%maximum_depth)) then
         allocate(integral)
      else
         allocate(bounded_integral)
         bounded_integral%minimum_depth = self%minimum_depth
         bounded_integral%maximum_depth = self%maximum_depth
         integral => bounded_integral
      end if
      integral%average = self%average
      call model%add_child(integral, trim(self%output_name) // '_calculator')
      call integral%request_coupling(integral%id_input, self%link)
      link => integral%id_output%link
   end function vertical_integral_initialize

   function interior_temporal_mean(input, period, resolution, missing_value) result(expression)
      class (type_dependency_id), intent(inout), target :: input
      real(rk),                   intent(in)            :: period, resolution
      real(rk), optional,         intent(in)            :: missing_value
      type (type_interior_temporal_mean)                :: expression

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::interior_temporal_mean', &
         'Input variable has not been registered yet.')

      write (expression%output_name,'(i0,a,a,a,i0,a)') int(period), '_s_mean_', trim(input%link%name), '_at_', int(resolution), '_s_resolution'
      expression%link => input%link
      expression%n = nint(period / resolution)
      expression%period = period
      expression%use_incomplete_result = .not. present(missing_value)
      if (present(missing_value)) expression%missing_value = missing_value
   end function

   function horizontal_temporal_mean(input, period, resolution, missing_value) result(expression)
      class (type_horizontal_dependency_id), intent(inout), target :: input
      real(rk),                              intent(in)            :: period, resolution
      real(rk), optional,                    intent(in)            :: missing_value
      type (type_horizontal_temporal_mean)                         :: expression

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::horizontal_temporal_mean', &
         'Input variable has not been registered yet.')

      write (expression%output_name,'(i0,a,a,a,i0,a)') int(period), '_s_mean_', trim(input%link%name), '_at_', int(resolution), '_s_resolution'
      expression%link => input%link
      expression%n = nint(period / resolution)
      expression%period = period
      expression%use_incomplete_result = .not. present(missing_value)
      if (present(missing_value)) expression%missing_value = missing_value
   end function

   function horizontal_temporal_maximum(input, period, resolution, missing_value) result(expression)
      class (type_horizontal_dependency_id), intent(inout), target :: input
      real(rk),                              intent(in)            :: period, resolution
      real(rk), optional,                    intent(in)            :: missing_value
      type (type_horizontal_temporal_maximum)                      :: expression

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::horizontal_temporal_max', &
         'Input variable has not been registered yet.')

      write (expression%output_name,'(i0,a,a,a,i0,a)') int(period), '_s_max_', trim(input%link%name), '_at_', int(resolution), '_s_resolution'
      expression%link => input%link
      expression%n = nint(period / resolution)
      expression%period = period
      if (present(missing_value)) expression%missing_value = missing_value
   end function

   function interior_temporal_mean_initialize(self, model) result(link)
      class (type_interior_temporal_mean), intent(inout) :: self
      class (type_base_model),             intent(inout) :: model
      type (type_link), pointer :: link

      integer :: ibin
      character(len=10) :: strindex

      allocate(self%history(self%n + 1))
      do ibin = 1, size(self%history)
         write (strindex,'(i0)') ibin
         call model%add_interior_variable(get_safe_name(trim(self%output_name)) // "_bin" // trim(strindex), link=self%history(ibin)%link, source=source_expression, output=output_none)
         self%history(ibin)%link%target%prefill_value = 0.0_rk
         self%history(ibin)%link%target%part_of_state = .true.
      end do
      call model%add_interior_variable(get_safe_name(trim(self%output_name)) // "_last", link=self%previous_value%link, source=source_expression, output=output_none)
      call model%add_interior_variable(get_safe_name(trim(self%output_name)) // "_last_exact_mean", link=self%last_exact_mean%link, source=source_expression, output=output_none)
      call model%add_interior_variable(get_safe_name(trim(self%output_name)) // "_mean", link=self%mean%link, source=source_expression, output=output_none)
      call model%add_scalar_variable(get_safe_name(trim(self%output_name)) // "_last_time", link=self%previous_time%link, source=source_expression, output=output_none)
      call model%add_scalar_variable(get_safe_name(trim(self%output_name)) // "_start_time", link=self%start_time%link, source=source_expression, output=output_none)
      call model%add_scalar_variable(get_safe_name(trim(self%output_name)) // "_icurrent", link=self%icurrent%link, source=source_expression, output=output_none)
      self%last_exact_mean%link%target%prefill_value = 0.0_rk
      self%mean%link%target%prefill_value = self%missing_value
      self%start_time%link%target%prefill_value = -huge(self%start_time%p)
      self%icurrent%link%target%prefill_value = 1.0_rk
      self%previous_value%link%target%part_of_state = .true.
      self%previous_time%link%target%part_of_state = .true.
      self%start_time%link%target%part_of_state = .true.
      self%icurrent%link%target%part_of_state = .true.
      self%last_exact_mean%link%target%part_of_state = .true.
      link => self%mean%link
   end function interior_temporal_mean_initialize

   subroutine interior_temporal_mean_set_data(self, interior_store, horizontal_store, scalar_store, seconds_per_time_unit)
      class (type_interior_temporal_mean), intent(inout) :: self
      real(rke), target _CONTIGUOUS_, dimension(_PREARG_LOCATION_DIMENSIONS_ 0:)            :: interior_store
      real(rke), target _CONTIGUOUS_, dimension(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ 0:) :: horizontal_store
      real(rke), target _CONTIGUOUS_, dimension(0:)                                         :: scalar_store
      real(rke), intent(in)                                                                 :: seconds_per_time_unit

      integer :: ibin

      self%in = self%link%target%catalog_index
      self%period = self%period / seconds_per_time_unit
      do ibin = 1, size(self%history)
         self%history(ibin)%p => interior_store(_PREARG_LOCATION_DIMENSIONS_ self%history(ibin)%link%target%store_index)
      end do
      self%previous_value%p => interior_store(_PREARG_LOCATION_DIMENSIONS_ self%previous_value%link%target%store_index)
      self%last_exact_mean%p => interior_store(_PREARG_LOCATION_DIMENSIONS_ self%last_exact_mean%link%target%store_index)
      self%mean%p => interior_store(_PREARG_LOCATION_DIMENSIONS_ self%mean%link%target%store_index)
      self%previous_time%p => scalar_store(self%previous_time%link%target%store_index)
      self%start_time%p => scalar_store(self%start_time%link%target%store_index)
      self%icurrent%p => scalar_store(self%icurrent%link%target%store_index)
   end subroutine

   subroutine interior_temporal_mean_update(self, time, value _POSTARG_LOCATION_RANGE_)
      class (type_interior_temporal_mean), intent(inout) :: self
      real(rke),                           intent(in)    :: time
      real(rke) _ATTRIBUTES_GLOBAL_,       intent(in)    :: value
      _DECLARE_ARGUMENTS_LOCATION_RANGE_

      real(rke) :: dt, w, dt_bin, scale, bin_end_time
      integer :: icurrent, icurrentbin, ioldest
      _DECLARE_LOCATION_

      ! Note that all array processing below uses explicit loops in order to respect
      ! any limits on the active domain given by the _LOCATION_RANGE_ argument.

      dt_bin = self%period / self%n

      if (self%start_time%p == -huge(self%start_time%p)) then
         self%start_time%p = time
         self%previous_time%p = time
      end if

      icurrent = self%icurrent%p
      icurrentbin = mod(icurrent - 1, self%n + 1) + 1

      do
         ioldest = mod(icurrent, self%n + 1) + 1
         bin_end_time = self%start_time%p + dt_bin * icurrent
         if (bin_end_time > time) exit

         ! Linearly interpolate to value at end-of-bin time and add that to the current bin
         dt = bin_end_time - self%previous_time%p
         w = dt / (time - self%previous_time%p)   ! weight for current time (leaving 1-w for previous time)
         scale = dt / self%period
         _BEGIN_GLOBAL_LOOP_
            self%history(icurrentbin)%p _INDEX_LOCATION_  = self%history(icurrentbin)%p _INDEX_LOCATION_  &
               + ((1._rke - 0.5_rke * w) * self%previous_value%p _INDEX_LOCATION_ + 0.5_rke * w * value _INDEX_LOCATION_) &
               * scale
         _END_GLOBAL_LOOP_

         if (icurrent > self%n) then
            ! We already had a complete history (bins covering the full window size). Add the newly full bin, subtract the oldest bin
            _BEGIN_GLOBAL_LOOP_
               self%last_exact_mean%p _INDEX_LOCATION_ = self%last_exact_mean%p _INDEX_LOCATION_ &
                  - self%history(ioldest)%p _INDEX_LOCATION_  + self%history(icurrentbin)%p _INDEX_LOCATION_ 
            _END_GLOBAL_LOOP_
         else
            ! History is incomplete - just add newly filled bin
            _BEGIN_GLOBAL_LOOP_
               self%last_exact_mean%p _INDEX_LOCATION_ = self%last_exact_mean%p _INDEX_LOCATION_  &
                  + self%history(icurrentbin)%p _INDEX_LOCATION_ 
            _END_GLOBAL_LOOP_
         end if

         ! Update previous time and value to match end of current bin
         self%previous_time%p = bin_end_time
         _BEGIN_GLOBAL_LOOP_
            self%previous_value%p _INDEX_LOCATION_ = (1._rke - w) * self%previous_value%p _INDEX_LOCATION_ + w * value _INDEX_LOCATION_
         _END_GLOBAL_LOOP_

         ! Move to next bin: update indices, end time and empty newly current bin
         icurrent = icurrent + 1
         icurrentbin = mod(icurrent - 1, self%n + 1) + 1
         self%icurrent%p = icurrent
         _BEGIN_GLOBAL_LOOP_
            self%history(icurrentbin)%p _INDEX_LOCATION_  = 0
         _END_GLOBAL_LOOP_
      end do

      ! Compute average of previous and current value, multiply by time difference, pre-divide by window size, and add to current bin.
      scale = 0.5_rke * (time - self%previous_time%p) / self%period
      _BEGIN_GLOBAL_LOOP_
         self%history(icurrentbin)%p _INDEX_LOCATION_  = self%history(icurrentbin)%p _INDEX_LOCATION_  &
            + scale * (self%previous_value%p _INDEX_LOCATION_ + value _INDEX_LOCATION_)
      _END_GLOBAL_LOOP_

      if (icurrent > self%n) then
         ! We have a full history covering at least one window size. Update the running mean.
         ! The result is an approximation that assumes linear change over the period covered by the oldest bin.
         scale = (time - bin_end_time + dt_bin) / dt_bin
         _BEGIN_GLOBAL_LOOP_
            self%mean%p _INDEX_LOCATION_ = self%last_exact_mean%p _INDEX_LOCATION_ &
               + self%history(icurrentbin)%p _INDEX_LOCATION_  &
               - scale * self%history(ioldest)%p _INDEX_LOCATION_ 
         _END_GLOBAL_LOOP_
      elseif (self%use_incomplete_result) then
         if (self%start_time%p == time) then
            ! No results just - just the current (first) point in time
            _BEGIN_GLOBAL_LOOP_
               self%mean%p _INDEX_LOCATION_ = value _INDEX_LOCATION_
            _END_GLOBAL_LOOP_
         else
            ! Use average so far. The integral has been pre-divided by self%period.
            ! Undo this and divide instead by time integrated so far
            scale = self%period / (time - self%start_time%p)
            _BEGIN_GLOBAL_LOOP_
               self%mean%p _INDEX_LOCATION_ = scale * (self%last_exact_mean%p _INDEX_LOCATION_ &
                  + self%history(icurrentbin)%p _INDEX_LOCATION_)
            _END_GLOBAL_LOOP_
         end if
      end if

      ! Store current time and value to enable linear interpolation in subsequent call.
      self%previous_time%p = time
      _BEGIN_GLOBAL_LOOP_
         self%previous_value%p _INDEX_LOCATION_ = value _INDEX_LOCATION_
      _END_GLOBAL_LOOP_
   end subroutine

   function horizontal_temporal_mean_initialize(self, model) result(link)
      class (type_horizontal_temporal_mean), intent(inout) :: self
      class (type_base_model),               intent(inout) :: model
      type (type_link), pointer :: link

      integer :: ibin
      character(len=10) :: strindex

      allocate(self%history(self%n + 1))
      do ibin = 1, size(self%history)
         write (strindex,'(i0)') ibin
         call model%add_horizontal_variable(get_safe_name(trim(self%output_name)) // "_bin" // trim(strindex), link=self%history(ibin)%link, source=source_expression, output=output_none)
         self%history(ibin)%link%target%prefill_value = 0.0_rk
         self%history(ibin)%link%target%part_of_state = .true.
      end do
      call model%add_horizontal_variable(get_safe_name(trim(self%output_name)) // "_last", link=self%previous_value%link, source=source_expression, output=output_none)
      call model%add_horizontal_variable(get_safe_name(trim(self%output_name)) // "_last_exact_mean", link=self%last_exact_mean%link, source=source_expression, output=output_none)
      call model%add_horizontal_variable(get_safe_name(trim(self%output_name)) // "_mean", link=self%mean%link, source=source_expression, output=output_none)
      call model%add_scalar_variable(get_safe_name(trim(self%output_name)) // "_last_time", link=self%previous_time%link, source=source_expression, output=output_none)
      call model%add_scalar_variable(get_safe_name(trim(self%output_name)) // "_start_time", link=self%start_time%link, source=source_expression, output=output_none)
      call model%add_scalar_variable(get_safe_name(trim(self%output_name)) // "_icurrent", link=self%icurrent%link, source=source_expression, output=output_none)
      self%last_exact_mean%link%target%prefill_value = 0.0_rk
      self%mean%link%target%prefill_value = self%missing_value
      self%start_time%link%target%prefill_value = -huge(self%start_time%p)
      self%icurrent%link%target%prefill_value = 1.0_rk
      self%previous_value%link%target%part_of_state = .true.
      self%previous_time%link%target%part_of_state = .true.
      self%start_time%link%target%part_of_state = .true.
      self%icurrent%link%target%part_of_state = .true.
      self%last_exact_mean%link%target%part_of_state = .true.
      link => self%mean%link
   end function horizontal_temporal_mean_initialize

   subroutine horizontal_temporal_mean_set_data(self, interior_store, horizontal_store, scalar_store, seconds_per_time_unit)
      class (type_horizontal_temporal_mean), intent(inout) :: self
      real(rke), target _CONTIGUOUS_, dimension(_PREARG_LOCATION_DIMENSIONS_ 0:)            :: interior_store
      real(rke), target _CONTIGUOUS_, dimension(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ 0:) :: horizontal_store
      real(rke), target _CONTIGUOUS_, dimension(0:)                                         :: scalar_store
      real(rke), intent(in)                                                                 :: seconds_per_time_unit

      integer :: ibin

      self%in = self%link%target%catalog_index
      self%period = self%period / seconds_per_time_unit
      do ibin = 1, size(self%history)
         self%history(ibin)%p => horizontal_store(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%history(ibin)%link%target%store_index)
      end do
      self%previous_value%p => horizontal_store(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%previous_value%link%target%store_index)
      self%last_exact_mean%p => horizontal_store(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%last_exact_mean%link%target%store_index)
      self%mean%p => horizontal_store(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%mean%link%target%store_index)
      self%previous_time%p => scalar_store(self%previous_time%link%target%store_index)
      self%start_time%p => scalar_store(self%start_time%link%target%store_index)
      self%icurrent%p => scalar_store(self%icurrent%link%target%store_index)
   end subroutine

   subroutine horizontal_temporal_mean_update(self, time, value _POSTARG_HORIZONTAL_LOCATION_RANGE_)
      class (type_horizontal_temporal_mean),    intent(inout) :: self
      real(rke),                                intent(in)    :: time
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, intent(in)    :: value
      _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_RANGE_

      real(rke) :: dt, w, dt_bin, scale, bin_end_time
      integer :: icurrent, icurrentbin, ioldest
      _DECLARE_HORIZONTAL_LOCATION_

      ! Note that all array processing below uses explicit loops in order to respect
      ! any limits on the active domain given by the _HORIZONTAL_LOCATION_RANGE_ argument.

      dt_bin = self%period / self%n

      if (self%start_time%p == -huge(self%start_time%p)) then
         self%start_time%p = time
         self%previous_time%p = time
      end if

      icurrent = self%icurrent%p
      icurrentbin = mod(icurrent - 1, self%n + 1) + 1

      do
         ioldest = mod(icurrent, self%n + 1) + 1
         bin_end_time = self%start_time%p + dt_bin * icurrent
         if (bin_end_time > time) exit

         ! Linearly interpolate to value at end-of-bin time and add that to the current bin
         dt = bin_end_time - self%previous_time%p
         w = dt / (time - self%previous_time%p)   ! weight for current time (leaving 1-w for previous time)
         scale = dt / self%period
         _BEGIN_OUTER_VERTICAL_LOOP_
            self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_  = self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_  &
               + ((1._rke - 0.5_rke * w) * self%previous_value%p _INDEX_HORIZONTAL_LOCATION_ + 0.5_rke * w * value _INDEX_HORIZONTAL_LOCATION_) &
               * scale
         _END_OUTER_VERTICAL_LOOP_

         if (icurrent > self%n) then
            ! We already had a complete history (bins covering the full window size). Add the newly full bin, subtract the oldest bin
            _BEGIN_OUTER_VERTICAL_LOOP_
               self%last_exact_mean%p _INDEX_HORIZONTAL_LOCATION_ = self%last_exact_mean%p _INDEX_HORIZONTAL_LOCATION_ &
                  - self%history(ioldest)%p _INDEX_HORIZONTAL_LOCATION_  + self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_ 
            _END_OUTER_VERTICAL_LOOP_
         else
            ! History is incomplete - just add newly filled bin
            _BEGIN_OUTER_VERTICAL_LOOP_
               self%last_exact_mean%p _INDEX_HORIZONTAL_LOCATION_ = self%last_exact_mean%p _INDEX_HORIZONTAL_LOCATION_  &
                  + self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_ 
            _END_OUTER_VERTICAL_LOOP_
         end if

         ! Update previous time and value to match end of current bin
         self%previous_time%p = bin_end_time
         _BEGIN_OUTER_VERTICAL_LOOP_
            self%previous_value%p _INDEX_HORIZONTAL_LOCATION_ = (1._rke - w) * self%previous_value%p _INDEX_HORIZONTAL_LOCATION_ + w * value _INDEX_HORIZONTAL_LOCATION_
         _END_OUTER_VERTICAL_LOOP_

         ! Move to next bin: update indices, end time and empty newly current bin
         icurrent = icurrent + 1
         icurrentbin = mod(icurrent - 1, self%n + 1) + 1
         self%icurrent%p = icurrent
         _BEGIN_OUTER_VERTICAL_LOOP_
            self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_ = 0
         _END_OUTER_VERTICAL_LOOP_
      end do

      ! Compute average of previous and current value, multiply by time difference, pre-divide by window size, and add to current bin.
      scale = 0.5_rke * (time - self%previous_time%p) / self%period
      _BEGIN_OUTER_VERTICAL_LOOP_
         self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_  = self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_  &
            + scale * (self%previous_value%p _INDEX_HORIZONTAL_LOCATION_ + value _INDEX_HORIZONTAL_LOCATION_)
      _END_OUTER_VERTICAL_LOOP_

      if (icurrent > self%n) then
         ! We have a full history covering at least one window size. Update the running mean.
         ! The result is an approximation that assumes linear change over the period covered by the oldest bin.
         scale = (time - bin_end_time + dt_bin) / dt_bin
         _BEGIN_OUTER_VERTICAL_LOOP_
            self%mean%p _INDEX_HORIZONTAL_LOCATION_ = self%last_exact_mean%p _INDEX_HORIZONTAL_LOCATION_ &
               + self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_  &
               - scale * self%history(ioldest)%p _INDEX_HORIZONTAL_LOCATION_
         _END_OUTER_VERTICAL_LOOP_
      elseif (self%use_incomplete_result) then
         if (self%start_time%p == time) then
            ! No results just - just the current (first) point in time
            _BEGIN_OUTER_VERTICAL_LOOP_
               self%mean%p _INDEX_HORIZONTAL_LOCATION_ = value _INDEX_HORIZONTAL_LOCATION_
            _END_OUTER_VERTICAL_LOOP_
         else
            ! Use average so far. The integral has been pre-divided by self%period.
            ! Undo this and divide instead by time integrated so far
            scale  = self%period / (time - self%start_time%p)
            _BEGIN_OUTER_VERTICAL_LOOP_
               self%mean%p _INDEX_HORIZONTAL_LOCATION_ = scale * (self%last_exact_mean%p _INDEX_HORIZONTAL_LOCATION_ &
                  + self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_)
            _END_OUTER_VERTICAL_LOOP_
         end if
      end if

      ! Store current time and value to enable linear interpolation in subsequent call.
      self%previous_time%p = time
      _BEGIN_OUTER_VERTICAL_LOOP_
         self%previous_value%p _INDEX_HORIZONTAL_LOCATION_ = value _INDEX_HORIZONTAL_LOCATION_
      _END_OUTER_VERTICAL_LOOP_
   end subroutine horizontal_temporal_mean_update

   function horizontal_temporal_maximum_initialize(self, model) result(link)
      class (type_horizontal_temporal_maximum), intent(inout) :: self
      class (type_base_model),                  intent(inout) :: model
      type (type_link), pointer :: link

      integer :: ibin
      character(len=10) :: strindex

      allocate(self%history(self%n))
      do ibin = 1, size(self%history)
         write (strindex,'(i0)') ibin
         call model%add_horizontal_variable(get_safe_name(trim(self%output_name)) // "_bin" // trim(strindex), link=self%history(ibin)%link, source=source_expression, output=output_none)
         self%history(ibin)%link%target%prefill_value = -huge(self%history(ibin)%link%target%prefill_value)
         self%history(ibin)%link%target%part_of_state = .true.
      end do
      call model%add_horizontal_variable(get_safe_name(trim(self%output_name)) // "_last", link=self%previous_value%link, source=source_expression, output=output_none)
      call model%add_horizontal_variable(get_safe_name(trim(self%output_name)) // "_max", link=self%maximum%link, source=source_expression, output=output_none)
      call model%add_scalar_variable(get_safe_name(trim(self%output_name)) // "_last_time", link=self%previous_time%link, source=source_expression, output=output_none)
      call model%add_scalar_variable(get_safe_name(trim(self%output_name)) // "_start_time", link=self%start_time%link, source=source_expression, output=output_none)
      call model%add_scalar_variable(get_safe_name(trim(self%output_name)) // "_icurrent", link=self%current%link, source=source_expression, output=output_none)
      self%maximum%link%target%prefill_value = self%missing_value
      self%start_time%link%target%prefill_value = -huge(self%start_time%p)
      self%current%link%target%prefill_value = 1.0_rk
      self%maximum%link%target%part_of_state = .true.
      self%previous_value%link%target%part_of_state = .true.
      self%previous_time%link%target%part_of_state = .true.
      self%start_time%link%target%part_of_state = .true.
      self%current%link%target%part_of_state = self%n > 1
      link => self%maximum%link
   end function horizontal_temporal_maximum_initialize

   subroutine horizontal_temporal_maximum_set_data(self, interior_store, horizontal_store, scalar_store, seconds_per_time_unit)
      class (type_horizontal_temporal_maximum), intent(inout) :: self
      real(rke), target _CONTIGUOUS_, dimension(_PREARG_LOCATION_DIMENSIONS_ 0:)            :: interior_store
      real(rke), target _CONTIGUOUS_, dimension(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ 0:) :: horizontal_store
      real(rke), target _CONTIGUOUS_, dimension(0:)                                         :: scalar_store
      real(rke), intent(in)                                                                 :: seconds_per_time_unit

      integer :: ibin

      self%in = self%link%target%catalog_index
      self%period = self%period / seconds_per_time_unit
      do ibin = 1, size(self%history)
         self%history(ibin)%p => horizontal_store(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%history(ibin)%link%target%store_index)
      end do
      self%previous_value%p =>horizontal_store(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%previous_value%link%target%store_index)
      self%maximum%p => horizontal_store(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%maximum%link%target%store_index)
      self%previous_time%p => scalar_store(self%previous_time%link%target%store_index)
      self%start_time%p => scalar_store(self%start_time%link%target%store_index)
      self%current%p => scalar_store(self%current%link%target%store_index)
   end subroutine

   subroutine horizontal_temporal_maximum_update(self, time, value _POSTARG_HORIZONTAL_LOCATION_RANGE_)
      class (type_horizontal_temporal_maximum), intent(inout) :: self
      real(rke),                                intent(in)    :: time
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, intent(in)    :: value
      _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_RANGE_

      integer :: ibin, icurrent, icurrentbin
      real(rke) :: w, bin_end_time
      _DECLARE_HORIZONTAL_LOCATION_

      ! Note that all array processing below uses explicit loops in order to respect
      ! any limits on the active domain given by the _HORIZONTAL_LOCATION_RANGE_ argument.

      if (self%start_time%p == -huge(self%start_time%p)) self%start_time%p = time

      icurrent = self%current%p
      icurrentbin = mod(icurrent - 1, self%n) + 1
      do
         bin_end_time = self%start_time%p + (self%period / self%n) * icurrent
         if (bin_end_time > time) exit

         ! Update previous time and value to match end of current bin. For the latter, linearly interpolate to value at end-of-bin time.
         w = (bin_end_time - self%previous_time%p) / (time - self%previous_time%p)   ! weight for current time (leaving 1-w for previous time)
         self%previous_time%p = bin_end_time
         _BEGIN_OUTER_VERTICAL_LOOP_
            self%previous_value%p _INDEX_HORIZONTAL_LOCATION_ = (1._rke - w) * self%previous_value%p _INDEX_HORIZONTAL_LOCATION_ + w * value _INDEX_HORIZONTAL_LOCATION_
         _END_OUTER_VERTICAL_LOOP_

         ! Complete the current bin by taking the maximum over its previous value and the value at end-of-bin time.
         _BEGIN_OUTER_VERTICAL_LOOP_
            self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_ = max(self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_, &
               self%previous_value%p _INDEX_HORIZONTAL_LOCATION_)
         _END_OUTER_VERTICAL_LOOP_

         if (icurrent >= self%n) then
            ! We have a complete history - compute the maximum over all bins
            _BEGIN_OUTER_VERTICAL_LOOP_
               self%maximum%p _INDEX_HORIZONTAL_LOCATION_ = self%history(1)%p _INDEX_HORIZONTAL_LOCATION_
            _END_OUTER_VERTICAL_LOOP_
            do ibin = 2, self%n
               _BEGIN_OUTER_VERTICAL_LOOP_
                  self%maximum%p _INDEX_HORIZONTAL_LOCATION_ = max(self%maximum%p _INDEX_HORIZONTAL_LOCATION_,  &
                     self%history(ibin)%p _INDEX_HORIZONTAL_LOCATION_)
               _END_OUTER_VERTICAL_LOOP_
            end do
         end if

         ! Move to next bin: update indices, end time and set maximum of newly current bin to current value (at start of bin)
         icurrent = icurrent + 1
         icurrentbin = mod(icurrent -1, self%n) + 1
         self%current%p = icurrent
         _BEGIN_OUTER_VERTICAL_LOOP_
            self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_ = self%previous_value%p _INDEX_HORIZONTAL_LOCATION_
         _END_OUTER_VERTICAL_LOOP_
      end do

      ! Update the maximum of the current bin
      _BEGIN_OUTER_VERTICAL_LOOP_
         self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_ = max(self%history(icurrentbin)%p _INDEX_HORIZONTAL_LOCATION_, &
            value _INDEX_HORIZONTAL_LOCATION_)
      _END_OUTER_VERTICAL_LOOP_

      ! Store current time and value to enable linear interpolation in subsequent call.
      self%previous_time%p = time
      _BEGIN_OUTER_VERTICAL_LOOP_
         self%previous_value%p _INDEX_HORIZONTAL_LOCATION_ = value _INDEX_HORIZONTAL_LOCATION_
      _END_OUTER_VERTICAL_LOOP_
   end subroutine

end module fabm_expressions

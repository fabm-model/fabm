#include "fabm_driver.h"
#include "fabm_private.h"

module fabm_builtin_time_filter
   use fabm_types
   use fabm_global_types

   implicit none

   private

   !public type_interior_temporal_mean, type_surface_temporal_mean, type_bottom_temporal_mean, type_surface_temporal_maximum
   public type_interior_temporal_mean, type_horizontal_temporal_mean, type_horizontal_temporal_maximum

   type, extends(type_global_model) :: type_interior_temporal_mean
      real(rk) :: window   ! Time period to average over (s)
      integer  :: n
      real(rk) :: missing_value = -2.e20_rk
      logical  :: use_incomplete_result = .false.

      type (type_global_interior_dependency_id) :: source
      type (type_global_interior_variable_id) :: mean
      type (type_global_interior_variable_id), private :: previous_value, last_exact_mean
      type (type_global_interior_variable_id), allocatable, private :: history(:)
      type (type_global_scalar_variable_id), private :: previous_time, start_time, icurrent
   contains
      procedure :: initialize => interior_temporal_mean_initialize
      procedure :: set_data   => interior_temporal_mean_set_data
      procedure :: update     => interior_temporal_mean_update
   end type

   type, extends(type_global_model) :: type_horizontal_temporal_mean
      real(rk) :: window   ! Time period to average over (s)
      integer  :: n
      real(rk) :: missing_value = -2.e20_rk
      logical  :: use_incomplete_result = .false.

      type (type_global_horizontal_dependency_id) :: source
      type (type_global_horizontal_variable_id) :: mean
      type (type_global_horizontal_variable_id), private :: previous_value, last_exact_mean
      type (type_global_horizontal_variable_id), allocatable, private :: history(:)
      type (type_global_scalar_variable_id), private :: previous_time, start_time, icurrent
   contains
      procedure :: initialize => horizontal_temporal_mean_initialize
      procedure :: set_data   => horizontal_temporal_mean_set_data
      procedure :: update     => horizontal_temporal_mean_update
   end type
   !
   !type, extends(type_global_model) :: type_bottom_temporal_mean
   !   type (type_bottom_diagnostic_variable_id) :: id_output
   !   type (type_horizontal_dependency_id)      :: id_input
   !contains
   !   procedure :: initialize => bottom_temporal_mean_initialize
   !   procedure :: set_data   => bottom_temporal_mean_set_data
   !   procedure :: update     => bottom_temporal_mean_update
   !end type

   type, extends(type_global_model) :: type_horizontal_temporal_maximum
      real(rk) :: window                     ! Time window to compute running maximum over (s)
      integer  :: n                          ! Number of bins to use to cover the period
      real(rk) :: missing_value = -2.e20_rk  ! Missing value to use until the simulation has covered the window size [period]

      type (type_global_horizontal_dependency_id) :: source
      type (type_global_horizontal_variable_id) :: maximum
      type (type_global_horizontal_variable_id), private :: previous_value
      type (type_global_horizontal_variable_id), allocatable, private :: history(:)
      type (type_global_scalar_variable_id), private :: previous_time, start_time, current
   contains
      procedure :: initialize => horizontal_temporal_maximum_initialize
      procedure :: set_data   => horizontal_temporal_maximum_set_data
      procedure :: update     => horizontal_temporal_maximum_update
   end type

contains

   subroutine interior_temporal_mean_initialize(self, configunit)
      class (type_interior_temporal_mean), intent(inout), target :: self
      integer,                             intent(in)            :: configunit

      integer :: ibin
      character(len=10) :: strindex

      if (self%user_created) then
         call self%get_parameter(self%window, 'window', 's', 'window size')
         call self%get_parameter(self%n, 'n', '', 'number of bins')
         call self%get_parameter(self%missing_value, 'missing_value', '', 'missing value to until the full window size has been covered', default=-2e20_rk)
      end if

      !call self%register_dependency(self%id_input, 'source', '', 'variable for which to compute running mean')
      !call self%register_diagnostic_variable(self%id_output, 'mean',  '', 'running mean')
      call self%add_interior_variable("source", link=self%source%link, presence=presence_external_required)

      allocate(self%history(self%n + 1))
      do ibin = 1, size(self%history)
         write (strindex,'(i0)') ibin
         call self%add_interior_variable("bin" // trim(strindex), link=self%history(ibin)%link, source=source_global, output=output_none)
         self%history(ibin)%link%target%prefill_value = 0.0_rk
         self%history(ibin)%link%target%part_of_state = .true.
      end do
      call self%add_interior_variable("last", link=self%previous_value%link, source=source_global, output=output_none)
      call self%add_interior_variable("last_exact_mean", link=self%last_exact_mean%link, source=source_global, output=output_none)
      call self%add_interior_variable("mean", link=self%mean%link, source=source_global, output=output_none)
      call self%add_scalar_variable("last_time", link=self%previous_time%link, source=source_global, output=output_none)
      call self%add_scalar_variable("start_time", link=self%start_time%link, source=source_global, output=output_none)
      call self%add_scalar_variable("icurrent", link=self%icurrent%link, source=source_global, output=output_none)
      self%last_exact_mean%link%target%prefill_value = 0.0_rk
      self%mean%link%target%prefill_value = self%missing_value
      self%start_time%link%target%prefill_value = -huge(self%start_time%p)
      self%icurrent%link%target%prefill_value = 1.0_rk
      self%previous_value%link%target%part_of_state = .true.
      self%previous_time%link%target%part_of_state = .true.
      self%start_time%link%target%part_of_state = .true.
      self%icurrent%link%target%part_of_state = .true.
      self%last_exact_mean%link%target%part_of_state = .true.
   end subroutine

   subroutine horizontal_temporal_mean_initialize(self, configunit)
      class (type_horizontal_temporal_mean), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      integer :: ibin
      character(len=10) :: strindex

      if (self%user_created) then
         call self%get_parameter(self%window, 'window', 's', 'window size')
         call self%get_parameter(self%n, 'n', '', 'number of bins')
         call self%get_parameter(self%missing_value, 'missing_value', '', 'missing value to use until the full window size has been covered', default=-2e20_rk)
      end if

      !call self%register_dependency(self%id_input, 'source', '', 'variable for which to compute running mean')
      !call self%register_diagnostic_variable(self%id_output, 'mean',  '', 'running mean')
      call self%add_horizontal_variable("source", link=self%source%link, presence=presence_external_required)

      allocate(self%history(self%n + 1))
      do ibin = 1, size(self%history)
         write (strindex,'(i0)') ibin
         call self%add_horizontal_variable("bin" // trim(strindex), link=self%history(ibin)%link, source=source_global, output=output_none)
         self%history(ibin)%link%target%prefill_value = 0.0_rk
         self%history(ibin)%link%target%part_of_state = .true.
      end do
      call self%add_horizontal_variable("last", link=self%previous_value%link, source=source_global, output=output_none)
      call self%add_horizontal_variable("last_exact_mean", link=self%last_exact_mean%link, source=source_global, output=output_none)
      call self%add_horizontal_variable("mean", link=self%mean%link, source=source_global, output=output_none)
      call self%add_scalar_variable("last_time", link=self%previous_time%link, source=source_global, output=output_none)
      call self%add_scalar_variable("start_time", link=self%start_time%link, source=source_global, output=output_none)
      call self%add_scalar_variable("icurrent", link=self%icurrent%link, source=source_global, output=output_none)
      self%last_exact_mean%link%target%prefill_value = 0.0_rk
      self%mean%link%target%prefill_value = self%missing_value
      self%start_time%link%target%prefill_value = -huge(self%start_time%p)
      self%icurrent%link%target%prefill_value = 1.0_rk
      self%previous_value%link%target%part_of_state = .true.
      self%previous_time%link%target%part_of_state = .true.
      self%start_time%link%target%part_of_state = .true.
      self%icurrent%link%target%part_of_state = .true.
      self%last_exact_mean%link%target%part_of_state = .true.
   end subroutine
   !
   !subroutine bottom_temporal_mean_initialize(self, configunit)
   !   class (type_bottom_temporal_mean), intent(inout), target :: self
   !   integer,                           intent(in)            :: configunit
   !
   !   real(rk) :: window, missing_value
   !   integer :: n
   !
   !   call self%get_parameter(self%window, 'window', 's', 'window size')
   !   call self%get_parameter(self%n, 'n', '', 'number of bins')
   !   call self%get_parameter(self%missing_value, 'missing_value', '', 'missing value to use until the full window size has been covered', default=-2e20_rk)
   !
   !   call self%register_dependency(self%id_input, 'source', '', 'variable for which to compute running mean')
   !   call self%register_diagnostic_variable(self%id_output, 'mean',  '', 'running mean')
   !end subroutine

   subroutine horizontal_temporal_maximum_initialize(self, configunit)
      class (type_horizontal_temporal_maximum), intent(inout), target :: self
      integer,                               intent(in)            :: configunit

      integer :: ibin
      character(len=10) :: strindex

      if (self%user_created) then
         call self%get_parameter(self%window, 'window', 's', 'window size')
         call self%get_parameter(self%n, 'n', '', 'number of bins')
         call self%get_parameter(self%missing_value, 'missing_value', '', 'missing value to use until the full window size has been covered', default=-2e20_rk)
      end if

      !call self%register_dependency(self%id_input, 'source', '', 'variable for which to compute running maximum')
      call self%add_horizontal_variable("source", link=self%source%link, presence=presence_external_required)

      allocate(self%history(self%n))
      do ibin = 1, size(self%history)
         write (strindex,'(i0)') ibin
         call self%add_horizontal_variable("bin" // trim(strindex), link=self%history(ibin)%link, source=source_global, output=output_none)
         self%history(ibin)%link%target%prefill_value = -huge(self%history(ibin)%link%target%prefill_value)
         self%history(ibin)%link%target%part_of_state = .true.
      end do
      call self%add_horizontal_variable("last", link=self%previous_value%link, source=source_global, output=output_none)
      call self%add_horizontal_variable("max", link=self%maximum%link, source=source_global, output=output_none)
      call self%add_scalar_variable("last_time", link=self%previous_time%link, source=source_global, output=output_none)
      call self%add_scalar_variable("start_time", link=self%start_time%link, source=source_global, output=output_none)
      call self%add_scalar_variable("icurrent", link=self%current%link, source=source_global, output=output_none)
      self%maximum%link%target%prefill_value = self%missing_value
      self%start_time%link%target%prefill_value = -huge(self%start_time%p)
      self%current%link%target%prefill_value = 1.0_rk
      self%maximum%link%target%part_of_state = .true.
      self%previous_value%link%target%part_of_state = .true.
      self%previous_time%link%target%part_of_state = .true.
      self%start_time%link%target%part_of_state = .true.
      self%current%link%target%part_of_state = self%n > 1
   end subroutine

   subroutine interior_temporal_mean_set_data(self, store, seconds_per_time_unit)
      class (type_interior_temporal_mean), intent(inout) :: self
      type (type_store), target                          :: store
      real(rke), optional, intent(in)                    :: seconds_per_time_unit

      integer :: ibin

      if (.not. present(seconds_per_time_unit)) call self%fatal_error('interior_temporal_mean_set_data', 'host did not provide time information.')
      self%source%icatalog = self%source%link%target%catalog_index
      self%window = self%window / seconds_per_time_unit
      do ibin = 1, size(self%history)
         self%history(ibin)%p => store%interior(_PREARG_LOCATION_DIMENSIONS_ self%history(ibin)%link%target%store_index)
      end do
      self%previous_value%p => store%interior(_PREARG_LOCATION_DIMENSIONS_ self%previous_value%link%target%store_index)
      self%last_exact_mean%p => store%interior(_PREARG_LOCATION_DIMENSIONS_ self%last_exact_mean%link%target%store_index)
      self%mean%p => store%interior(_PREARG_LOCATION_DIMENSIONS_ self%mean%link%target%store_index)
      self%previous_time%p => store%scalar(self%previous_time%link%target%store_index)
      self%start_time%p => store%scalar(self%start_time%link%target%store_index)
      self%icurrent%p => store%scalar(self%icurrent%link%target%store_index)
   end subroutine

   subroutine interior_temporal_mean_update(self, catalog _POSTARG_LOCATION_RANGE_, time)
      class (type_interior_temporal_mean), intent(in) :: self
      type (type_catalog),                 intent(in) :: catalog
      _DECLARE_ARGUMENTS_LOCATION_RANGE_
      real(rke), optional,                 intent(in) :: time

      real(rke) _ATTRIBUTES_GLOBAL_, pointer :: value
      real(rke) :: dt, w, dt_bin, scale, bin_end_time
      integer :: icurrent, icurrentbin, ioldest
      _DECLARE_LOCATION_

      if (.not. present(time)) call self%fatal_error('interior_temporal_mean_update', 'host did not provide time information.')

      value => catalog%interior(self%source%icatalog)%p
      if (.not. associated(value)) call self%fatal_error('interior_temporal_mean_update', 'source pointer not associated.')

      ! Note that all array processing below uses explicit loops in order to respect
      ! any limits on the active domain given by the _LOCATION_RANGE_ argument.

      dt_bin = self%window / self%n

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
         scale = dt / self%window
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
      scale = 0.5_rke * (time - self%previous_time%p) / self%window
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
            ! Use average so far. The integral has been pre-divided by self%window.
            ! Undo this and divide instead by time integrated so far
            scale = self%window / (time - self%start_time%p)
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

   subroutine horizontal_temporal_mean_set_data(self, store, seconds_per_time_unit)
      class (type_horizontal_temporal_mean), intent(inout) :: self
      type (type_store), target                            :: store
      real(rke), optional, intent(in)                      :: seconds_per_time_unit

      integer :: ibin

      if (.not. present(seconds_per_time_unit)) call self%fatal_error('horizontal_temporal_mean_set_data', 'host did not provide time information.')
      self%source%icatalog = self%source%link%target%catalog_index
      self%window = self%window / seconds_per_time_unit
      do ibin = 1, size(self%history)
         self%history(ibin)%p => store%horizontal(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%history(ibin)%link%target%store_index)
      end do
      self%previous_value%p => store%horizontal(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%previous_value%link%target%store_index)
      self%last_exact_mean%p => store%horizontal(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%last_exact_mean%link%target%store_index)
      self%mean%p => store%horizontal(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%mean%link%target%store_index)
      self%previous_time%p => store%scalar(self%previous_time%link%target%store_index)
      self%start_time%p => store%scalar(self%start_time%link%target%store_index)
      self%icurrent%p => store%scalar(self%icurrent%link%target%store_index)
   end subroutine

   subroutine horizontal_temporal_mean_update(self, catalog _POSTARG_LOCATION_RANGE_, time)
      class (type_horizontal_temporal_mean), intent(in) :: self
      type (type_catalog),                   intent(in) :: catalog
      _DECLARE_ARGUMENTS_LOCATION_RANGE_
      real(rke), optional,                   intent(in) :: time

      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: value
      real(rke) :: dt, w, dt_bin, scale, bin_end_time
      integer :: icurrent, icurrentbin, ioldest
      _DECLARE_HORIZONTAL_LOCATION_

      if (.not. present(time)) call self%fatal_error('horizontal_temporal_mean_update', 'host did not provide time information.')

      value => catalog%horizontal(self%source%icatalog)%p
      if (.not. associated(value)) call self%fatal_error('horizontal_temporal_mean_update', 'source pointer not associated.')

      ! Note that all array processing below uses explicit loops in order to respect
      ! any limits on the active domain given by the _LOCATION_RANGE_ argument.

      dt_bin = self%window / self%n

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
         scale = dt / self%window
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
      scale = 0.5_rke * (time - self%previous_time%p) / self%window
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
            ! Use average so far. The integral has been pre-divided by self%window.
            ! Undo this and divide instead by time integrated so far
            scale  = self%window / (time - self%start_time%p)
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

   subroutine horizontal_temporal_maximum_set_data(self, store, seconds_per_time_unit)
      class (type_horizontal_temporal_maximum), intent(inout) :: self
      type (type_store), target                               :: store
      real(rke), optional, intent(in)                         :: seconds_per_time_unit

      integer :: ibin

      if (.not. present(seconds_per_time_unit)) call self%fatal_error('horizontal_temporal_maximum_set_data', 'host did not provide time information.')
      self%source%icatalog = self%source%link%target%catalog_index
      self%window = self%window / seconds_per_time_unit
      do ibin = 1, size(self%history)
         self%history(ibin)%p => store%horizontal(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%history(ibin)%link%target%store_index)
      end do
      self%previous_value%p =>store%horizontal(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%previous_value%link%target%store_index)
      self%maximum%p => store%horizontal(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%maximum%link%target%store_index)
      self%previous_time%p => store%scalar(self%previous_time%link%target%store_index)
      self%start_time%p => store%scalar(self%start_time%link%target%store_index)
      self%current%p => store%scalar(self%current%link%target%store_index)
   end subroutine

   subroutine horizontal_temporal_maximum_update(self, catalog _POSTARG_LOCATION_RANGE_, time)
      class (type_horizontal_temporal_maximum), intent(in) :: self
      type (type_catalog),                      intent(in) :: catalog
      _DECLARE_ARGUMENTS_LOCATION_RANGE_
      real(rke), optional,                      intent(in) :: time

      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: value
      integer :: ibin, icurrent, icurrentbin
      real(rke) :: w, bin_end_time
      _DECLARE_HORIZONTAL_LOCATION_

      if (.not. present(time)) call self%fatal_error('horizontal_temporal_maximum_update', 'host did not provide time information.')

      value => catalog%horizontal(self%source%icatalog)%p
      if (.not. associated(value)) call self%fatal_error('horizontal_temporal_maximum_update', 'source pointer not associated.')

      ! Note that all array processing below uses explicit loops in order to respect
      ! any limits on the active domain given by the _LOCATION_RANGE_ argument.

      if (self%start_time%p == -huge(self%start_time%p)) self%start_time%p = time

      icurrent = self%current%p
      icurrentbin = mod(icurrent - 1, self%n) + 1
      do
         bin_end_time = self%start_time%p + (self%window / self%n) * icurrent
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

end module

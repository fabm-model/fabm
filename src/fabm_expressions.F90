#include "fabm_driver.h"

! This module define standard expressions that can be used by biogeochemical models.

module fabm_expressions

   use fabm_types
   use fabm_driver
   use fabm_builtin_depth_integral
   use fabm_builtin_time_filter

   implicit none

   private

   public temporal_mean, temporal_maximum, vertical_mean, vertical_integral

   type, extends(type_interior_expression) :: type_interior_temporal_mean_expression
      real(rk) :: period   ! Time period to average over (s)
      integer  :: n
      real(rk) :: missing_value = -2.e20_rk
      logical  :: use_incomplete_result = .false.

      type (type_link), pointer :: source => null()
   contains
      procedure :: resolve => interior_temporal_mean_resolve
   end type

   type, extends(type_horizontal_expression) :: type_horizontal_temporal_mean_expression
      real(rk) :: period   ! Time period to average over (s)
      integer  :: n
      real(rk) :: missing_value = 0.0_rk
      logical  :: use_incomplete_result = .false.

      type (type_link), pointer :: source => null()
   contains
      procedure :: resolve => horizontal_temporal_mean_resolve
   end type

   type, extends(type_horizontal_expression) :: type_horizontal_temporal_maximum_expression
      real(rk) :: period                     ! Time window to compute running maximum over (s)
      integer  :: n                          ! Number of bins to use to cover the period
      real(rk) :: missing_value = -2.e20_rk  ! Missing value to use until the simulation has covered the window size [period]

      type (type_link), pointer :: source => null()
   contains
      procedure :: resolve => horizontal_temporal_maximum_resolve
   end type

   type, extends(type_horizontal_expression) :: type_vertical_integral
      real(rk) :: minimum_depth = 0.0_rk        ! Depth below surface in m (positive)
      real(rk) :: maximum_depth = huge(1.0_rk)  ! Depth below surface in m (positive)
      logical  :: average       = .false.       ! Whether to divide the depth integral by water depth, thus computing the vertical average

      type (type_link), pointer :: source => null()
   contains
      procedure :: resolve => vertical_integral_resolve
   end type

   interface temporal_mean
      module procedure interior_temporal_mean
      module procedure horizontal_temporal_mean
   end interface

   interface temporal_maximum
      module procedure horizontal_temporal_maximum
   end interface

contains

   function vertical_integral(input, minimum_depth, maximum_depth, average) result(base_expression)
      class (type_dependency_id), intent(inout), target :: input
      real(rk), optional,         intent(in)            :: minimum_depth, maximum_depth
      logical,  optional,         intent(in)            :: average
      class (type_horizontal_expression), pointer       :: base_expression

      character(len=attribute_length)        :: postfix
      type (type_vertical_integral), pointer :: expression

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::vertical_integral', &
         'Input variable has not been registered yet.')

      ! Create a name for the expression
      if (present(minimum_depth) .and. present(maximum_depth)) then
         if (minimum_depth > maximum_depth) call fatal_error('fabm_expressions::vertical_integral', &
            'Minimum depth exceeds maximum depth.')
         write (postfix,'(a,i0,a,i0,a)') '_between_', int(minimum_depth), '_m_and_', int(maximum_depth), '_m'
      elseif (present(minimum_depth)) then
         write (postfix,'(a,i0,a)') '_below_', int(minimum_depth), '_m'
      elseif (present(maximum_depth)) then
         write (postfix,'(a,i0,a)') '_above_', int(maximum_depth), '_m'
      else
         postfix = ''
      end if

      allocate(expression)
      if (present(average)) expression%average = average
      if (expression%average) then
         expression%output_name = 'vertical_mean_' // trim(input%link%name) // trim(postfix)
      else
         expression%output_name = 'integral_of_' // trim(input%link%name) // '_wrt_depth' // trim(postfix)
      end if
      expression%source => input%link
      if (present(minimum_depth)) expression%minimum_depth = minimum_depth
      if (present(maximum_depth)) expression%maximum_depth = maximum_depth
      base_expression => expression
   end function

   function vertical_mean(input, minimum_depth, maximum_depth) result(base_expression)
      class (type_dependency_id), intent(inout), target :: input
      real(rk), optional,        intent(in)             :: minimum_depth,maximum_depth
      class (type_horizontal_expression), pointer       :: base_expression
      base_expression => vertical_integral(input, minimum_depth, maximum_depth, average=.true.)
   end function

   function vertical_integral_resolve(self) result(link)
      class (type_vertical_integral), intent(inout) :: self
      type (type_link), pointer                     :: link

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
      call self%link%target%owner%add_child(integral, trim(self%output_name) // '_calculator')
      call integral%request_coupling(integral%id_input, self%source)
      integral%id_output%link%target%output = output_none
      link => integral%id_output%link
   end function

   function interior_temporal_mean(input, period, resolution, missing_value) result(base_expression)
      class (type_dependency_id), intent(inout), target :: input
      real(rk),                   intent(in)            :: period, resolution
      real(rk), optional,         intent(in)            :: missing_value
      class (type_interior_expression), pointer         :: base_expression

      type (type_interior_temporal_mean_expression), pointer :: expression

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::interior_temporal_mean', &
         'Input variable has not been registered yet.')

      allocate(expression)
      write (expression%output_name,'(i0,a,a,a,i0,a)') int(period), '_s_mean_', trim(input%link%name), '_at_', int(resolution), '_s_resolution'
      expression%source => input%link
      expression%n = nint(period / resolution)
      expression%period = period
      expression%use_incomplete_result = .not. present(missing_value)
      if (present(missing_value)) expression%missing_value = missing_value
      base_expression => expression
   end function

   function interior_temporal_mean_resolve(self) result(link)
      class (type_interior_temporal_mean_expression), intent(inout) :: self
      type (type_link), pointer                                     :: link

      class (type_interior_temporal_mean), pointer :: calculator

      allocate(calculator)
      calculator%window = self%period
      calculator%n = self%n
      calculator%missing_value = self%missing_value
      calculator%use_incomplete_result = self%use_incomplete_result
      call self%link%target%owner%add_child(calculator, trim(self%output_name) // '_calculator')
      call calculator%request_coupling(calculator%source%link, self%source)
      calculator%mean%link%target%output = output_none
      link => calculator%mean%link
   end function

   function horizontal_temporal_mean(input, period, resolution, missing_value) result(base_expression)
      class (type_horizontal_dependency_id), intent(inout), target :: input
      real(rk),                              intent(in)            :: period, resolution
      real(rk), optional,                    intent(in)            :: missing_value
      class (type_horizontal_expression), pointer                  :: base_expression

      type (type_horizontal_temporal_mean_expression), pointer :: expression

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::horizontal_temporal_mean', &
         'Input variable has not been registered yet.')

      allocate(expression)
      write (expression%output_name,'(i0,a,a,a,i0,a)') int(period), '_s_mean_', trim(input%link%name), '_at_', int(resolution), '_s_resolution'
      expression%source => input%link
      expression%n = nint(period / resolution)
      expression%period = period
      expression%use_incomplete_result = .not. present(missing_value)
      if (present(missing_value)) expression%missing_value = missing_value
      base_expression => expression
   end function

   function horizontal_temporal_mean_resolve(self) result(link)
      class (type_horizontal_temporal_mean_expression), intent(inout) :: self
      type (type_link), pointer                                       :: link

      class (type_horizontal_temporal_mean), pointer :: calculator

      allocate(calculator)
      calculator%window = self%period
      calculator%n = self%n
      calculator%missing_value = self%missing_value
      calculator%use_incomplete_result = self%use_incomplete_result
      call self%link%target%owner%add_child(calculator, trim(self%output_name) // '_calculator')
      call calculator%request_coupling(calculator%source%link, self%source)
      calculator%mean%link%target%output = output_none
      link => calculator%mean%link
   end function

   function horizontal_temporal_maximum(input, period, resolution, missing_value) result(base_expression)
      class (type_horizontal_dependency_id), intent(inout), target :: input
      real(rk),                              intent(in)            :: period, resolution
      real(rk), optional,                    intent(in)            :: missing_value
      class (type_horizontal_expression), pointer                  :: base_expression

      type (type_horizontal_temporal_maximum_expression), pointer :: expression

      if (.not. associated(input%link)) call fatal_error('fabm_expressions::horizontal_temporal_max', &
         'Input variable has not been registered yet.')

      allocate(expression)
      write (expression%output_name,'(i0,a,a,a,i0,a)') int(period), '_s_max_', trim(input%link%name), '_at_', int(resolution), '_s_resolution'
      expression%source => input%link
      expression%n = nint(period / resolution)
      expression%period = period
      if (present(missing_value)) expression%missing_value = missing_value
      base_expression => expression
   end function

   function horizontal_temporal_maximum_resolve(self) result(link)
      class (type_horizontal_temporal_maximum_expression), intent(inout) :: self
      type (type_link), pointer                                          :: link

      class (type_horizontal_temporal_maximum), pointer :: calculator

      allocate(calculator)
      calculator%window = self%period
      calculator%n = self%n
      calculator%missing_value = self%missing_value
      call self%link%target%owner%add_child(calculator, trim(self%output_name) // '_calculator')
      call calculator%request_coupling(calculator%source%link, self%source)
      calculator%maximum%link%target%output = output_none
      link => calculator%maximum%link
   end function

end module fabm_expressions

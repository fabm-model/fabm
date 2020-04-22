#include "fabm_driver.h"
   
module nonlocal

   use fabm_types

   implicit none

   private 

   ! type_vertical_distribution: module that computes a "vertical distribution" defined in terms
   ! of weights applied to different depth levels when computing a vertical integral (e.g., of prey).
   ! In this simple model, the vertical distribution is defined by an absolute upper depth limit (z_top),
   ! an absolute lower depth limit (z_bot), a weight to apply at the top of the domain (w_top) and a weight to
   ! apply at the bottom of the domain (w_bot). Linear interpolation is used to determine weights in between
   ! these levels. In the special case where both weights equal 1, a vertical integral computed using
   ! this presence distribution will in fact be the real vertical integral between z_top and z_bot.
   type, extends(type_base_model), public :: type_vertical_distribution
      type (type_diagnostic_variable_id) :: id_weights
      type (type_dependency_id)          :: id_depth

      real(rk) :: z_top, z_bot, w_top, w_bot
   contains
      procedure :: initialize => vertical_distribution_initialize
      procedure :: do         => vertical_distribution_do
   end type

   ! type_depth_integral: module that computes weighted depth integral of a specified variable.
   ! Weights are taken from an externally computed variable [id_target below],
   ! which typically is a diagnostic computed by a model of type type_vertical_distribution.
   ! The variable to be integrated is a dependency [id_target below] that must be coupled at run time.
   ! The resultant integral is a diagnostic that acts like a state variable, that is to say,
   ! other models can use it as if it were a state variable, and provide sources and sinks,
   ! and this model will then automagically distribute those sinks and sources again over their
   ! original depth-explicit source variable, using the appropriate weights.
   type, extends(type_base_model), public :: type_depth_integral
      type (type_bottom_diagnostic_variable_id) :: id_integral
      type (type_state_variable_id)             :: id_target
      type (type_dependency_id)                 :: id_weights
      type (type_dependency_id)                 :: id_thickness
   contains
      procedure :: initialize => depth_integral_initialize
      procedure :: do_column  => depth_integral_do_column
   end type

   type, extends(type_base_model) :: type_depth_integral_rate_distributor
      type (type_bottom_dependency_id) :: id_integral   ! Depth-integrated target variable
      type (type_bottom_dependency_id) :: id_sms        ! Depth-integrated sources-sinks of target variable
      type (type_state_variable_id)    :: id_target     ! Depth-explicit variable that should absorp the sources-sinks
      type (type_dependency_id)        :: id_weights    ! Weights for the vertical distribution of the sinks and sources
   contains
      procedure :: initialize => depth_integral_rate_distributor_initialize
      procedure :: do         => depth_integral_rate_distributor_do
   end type
   
contains

   subroutine vertical_distribution_initialize(self, configunit)
      class (type_vertical_distribution), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      call self%get_parameter(self%z_top,'z_top','m','upper limit of vertical distribution')
      call self%get_parameter(self%z_bot,'z_bot','m','lower limit of vertical distribution')
      call self%get_parameter(self%w_top,'w_top','-','weight at upper limit of vertical distribution')
      call self%get_parameter(self%w_bot,'w_bot','-','weight at lower limit of vertical distribution')

      call self%register_diagnostic_variable(self%id_weights,'weights','-','weights')
      call self%register_dependency(self%id_depth,standard_variables%depth)
   end subroutine vertical_distribution_initialize

   subroutine vertical_distribution_do(self, _ARGUMENTS_DO_)
      class (type_vertical_distribution), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: z, slope

      _LOOP_BEGIN_
         _GET_(self%id_depth,z)
         if (z < self%z_top .or. z > self%z_bot) then
            ! Outside specified domain - relative presence is zero.
            _SET_DIAGNOSTIC_(self%id_weights,0.0_rk)
         else
            ! Inside specified domain: linear interpolation between presence coefficients
            slope = (self%w_bot-self%w_top)/(self%z_bot-self%z_top)
            _SET_DIAGNOSTIC_(self%id_weights,self%w_top+(z-self%z_top)*slope)

            ! Here the weights could be adjusted for those cells that are only partially within the specified domain.
         end if
      _LOOP_END_
   end subroutine vertical_distribution_do

   subroutine depth_integral_initialize(self, configunit)
      class (type_depth_integral), intent(inout), target :: self
      integer,                     intent(in)            :: configunit

      class (type_depth_integral_rate_distributor), pointer :: rate_distributor

      call self%register_state_dependency(self%id_target, 'target', '', 'variable to depth-integrate')
      call self%register_dependency(self%id_weights,'weights','-','weights for vertical integration')
      call self%register_diagnostic_variable(self%id_integral,'result','','result',missing_value=0.0_rk, &
         act_as_state_variable=.true.,source=source_do_column)
      call self%register_dependency(self%id_thickness,standard_variables%cell_thickness)

      ! Create a child model that receives the intended rate of the change of the depth-integrated variable,
      ! and redistributes that rate of change over the pelagic variable that the integral was originally computed from.
      ! All dependencies of the child model can be resolved by coupling to our own variables.
      allocate(rate_distributor)
      call self%add_child(rate_distributor,'rate_distributor',configunit=-1)
      call rate_distributor%request_coupling(rate_distributor%id_target,'target')
      call rate_distributor%request_coupling(rate_distributor%id_weights,'weights')
      call rate_distributor%request_coupling(rate_distributor%id_integral,'result')
      call rate_distributor%request_coupling(rate_distributor%id_sms,'result_sms_tot')
   end subroutine depth_integral_initialize

   subroutine depth_integral_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_depth_integral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: local, weight, thickness, integral

      integral = 0.0_rk
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_target,local)
         _GET_(self%id_weights,weight)
         _GET_(self%id_thickness,thickness)
         integral = integral + local*weight*thickness
      _VERTICAL_LOOP_END_
      _SET_BOTTOM_DIAGNOSTIC_(self%id_integral,integral)
   end subroutine depth_integral_do_column

   subroutine depth_integral_rate_distributor_initialize(self, configunit)
      class (type_depth_integral_rate_distributor), intent(inout), target :: self
      integer,                                      intent(in)            :: configunit

      call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
      call self%register_dependency(self%id_weights,'weights','-','weights for vertical distribution')
      call self%register_dependency(self%id_integral,'integral','','depth-integrated target variable')
      call self%register_dependency(self%id_sms,'sms','','depth-integrated sources-sinks')
   end subroutine depth_integral_rate_distributor_initialize

   subroutine depth_integral_rate_distributor_do(self, _ARGUMENTS_DO_)
      class (type_depth_integral_rate_distributor), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: integral, integrated_sms, relative_change
      real(rk) :: local, weight

      _LOOP_BEGIN_
         ! First compute relative rate of change of depth-integrated target variable.
         _GET_BOTTOM_(self%id_integral,integral)
         _GET_BOTTOM_(self%id_sms,integrated_sms)
         if (integral /= 0.0_rk) then
            relative_change = integrated_sms / integral
         else
            relative_change = 0.0_rk
         end if

         ! Now distribute the same relative change across the column.
         _GET_(self%id_target,local)
         _GET_(self%id_weights,weight)
         _ADD_SOURCE_(self%id_target,relative_change * local * weight)
      _LOOP_END_
   end subroutine depth_integral_rate_distributor_do

end module
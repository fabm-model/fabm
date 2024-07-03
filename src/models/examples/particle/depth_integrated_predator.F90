#include "fabm_driver.h"

module examples_particle_depth_integrated_predator

   ! This module describes a predator that feeds across (part of) the water column.
   ! It is represented by a depth-integrated biomass and a prescribed vertical distribution.
   ! This would be appropriate for predators that move very fast in the vertical.
   !
   ! The predator has constant C:N:P stoichiometry
   ! Its elemental ratios are defined by parameters NC and PC below.
   ! It accepts one prey type that may have variable C:N:P stochiometry.
   ! At each point in time, ingested fluxes of C, N and P are calculated.
   ! Based on the most limiting of these, a growth rate is calculated.
   ! Unused ingested fluxes and dead biomass are sent to a coupled waste pool.
   !
   ! Extensions:
   ! * Different vertical distributions, including time-varying ones:
   !   create a new depth distribution type (your equivalent of type_vertical_depth_range
   !   defined in the fabm_builtin_depth_mapping module) and calculate your custom
   !   vertical distribution weights in its "do" routine. These can depend on any
   !   environmental input or parameter, e.g. temperature, light, time of day,
   !   depth, prey availability.
   ! * Multiple prey types: declare id_prey_c and the like as allocatable arrays,
   !   get number of prey types as parameter, allocate identifier arrays (id_prey_c etc.)
   !   and move all prey handling to the inside of a loop over all prey
   ! * Additional prey constituents (e.g., silicon, calcium carbonate): as for
   !   current treatment of N and P. If not incorporated in predator biomass,
   !   their ingested fluxes can be sent directly to the waste pool.
   ! * More complex dynamics for predation, predator population growth, etc.:
   !   modify logic in do_surface accordingly. Additional parameters should be added
   !   to the type_depth_integrated_predator type and their value retrieved from
   !   initialize. Additional environental inputs can be handled just like
   !   temperature. Any additional waste pools (e.g., to distinguish dissolved
   !   and particulate wastes) can be implemented analogous to the current waste.

   use fabm_types
   use fabm_particle
   use fabm_builtin_depth_mapping

   implicit none

   private

   type, extends(type_depth_integrated_particle), public :: type_depth_integrated_predator
      type (type_surface_state_variable_id)      :: id_c
      type (type_surface_dependency_id)          :: id_prey_c, id_prey_n, id_prey_p
      type (type_surface_state_variable_id)      :: id_waste_c, id_waste_n, id_waste_p
      type (type_surface_dependency_id)          :: id_temp
      type (type_surface_diagnostic_variable_id) :: id_net_growth, id_prey_loss_rate
      type (type_model_id)                       :: id_prey_int

      real(rk) :: clearance_rate
      real(rk) :: mortality
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

   ! Redfieldian N:C and P:C ratios of predator biomass
   real(rk), parameter :: NC = 16.0_rk/ 106.0_rk
   real(rk), parameter :: PC = 1.0_rk / 106.0_rk

contains

   subroutine initialize(self, configunit)
      class (type_depth_integrated_predator), intent(inout), target :: self
      integer,                                intent(in)            :: configunit

      class (type_vertical_depth_range), pointer :: depth_distribution

      ! This adds a dependency ("w") on the weights that define the vertical habitat
      ! Additionally, it adds a dependency on the vertical integral of the weights (id_w%integral),
      ! which can be used to convert between depth integrals and depth-averages
      call self%type_depth_integrated_particle%initialize(configunit)

      ! For now, specify a prescribed depth distribution
      allocate(depth_distribution)
      call self%add_child(depth_distribution, 'habitat')

      ! Link our own weights that specify the vertical habitat (variable "w", inherited
      ! from type_depth_integrated_particle) to those computed by the "habitat" child model
      call self%request_coupling('w', 'habitat/w')

      ! Predator biomass and its contribution to different elemental pools
      call self%register_state_variable(self%id_c, 'c', 'mmol C m-2', 'density')
      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=NC)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=PC)

      ! Parameters
      ! NB rates are in d-1 in fabm.yaml and scaled here to s-1 using the scale_factor argument
      call self%get_parameter(self%clearance_rate, 'clearance_rate', 'm3 d-1 mmol-1', 'clearance rate', scale_factor=1.0_rk / 86400.0_rk)
      call self%get_parameter(self%mortality, 'mortality', 'd-1', 'mortality', scale_factor=1.0_rk / 86400.0_rk)

      ! Diagnostics
      call self%register_diagnostic_variable(self%id_net_growth, 'net_growth', 'mmol C m-2 d-1', 'net population growth rate')
      call self%register_diagnostic_variable(self%id_prey_loss_rate, 'prey_loss_rate', 'd-1', 'specific prey loss rate')

      ! Depth-averaged dependencies
      call self%register_dependency(self%id_temp, 'temp', 'degrees_Celsius', 'depth-averaged temperature')
      call self%register_dependency(self%id_prey_c, 'prey_c', 'mmol C m-3', 'depth-averaged prey carbon')
      call self%register_dependency(self%id_prey_n, 'prey_n', 'mmol N m-3', 'depth-averaged prey nitrogen')
      call self%register_dependency(self%id_prey_p, 'prey_p', 'mmol P m-3', 'depth-averaged prey phosphorus')
      call self%register_state_dependency(self%id_waste_c, 'waste_c', 'mmol C m-2', 'depth-integrated carbon waste')
      call self%register_state_dependency(self%id_waste_n, 'waste_n', 'mmol N m-2', 'depth-integrated nitrogen waste')
      call self%register_state_dependency(self%id_waste_p, 'waste_p', 'mmol P m-2', 'depth-integrated phosphorus waste')

      ! Derive depth-averaged dependencies from depth-explicit sources
      ! * Environmental variables are typically depth-averaged over the predator habitat,
      !   as temperature is below (note average=.true.)
      ! * Prey concentrations are here depth-averaged as well; they will be multiplied with a clearance rate
      !   (volume searched per unit time per predator) to obtain ingestion rates.
      !   Prey destruction in handled separately by coupling to all prey state variables at once
      !   (see the call to register_mapped_model_dependency below)
      ! * Waste pools ar depth-integrated as they will receive depth-integrated fluxes of waste
      !   produced by the predator population. These fluxes will be vertically distributed
      !   accordingly to the predator's vertical distribution, i.e., the waste flux injected
      !   locally will be proportional to the local weight of the predators's
      !   vertical distribution.
      call self%request_mapped_coupling(self%id_temp, standard_variables%temperature, average=.true.)
      call self%request_mapped_coupling_to_model(self%id_prey_c, 'prey', standard_variables%total_carbon, average=.true.)
      call self%request_mapped_coupling_to_model(self%id_prey_n, 'prey', standard_variables%total_nitrogen, average=.true.)
      call self%request_mapped_coupling_to_model(self%id_prey_p, 'prey', standard_variables%total_phosphorus, average=.true.)
      call self%request_mapped_coupling_to_model(self%id_waste_c, 'waste', standard_variables%total_carbon)
      call self%request_mapped_coupling_to_model(self%id_waste_n, 'waste', standard_variables%total_nitrogen)
      call self%request_mapped_coupling_to_model(self%id_waste_p, 'waste', standard_variables%total_phosphorus)

      ! Access depth-integrated prey state that we will apply specific loss rates to.
      ! The local (depth-explicit) prey loss is expected to be proportional to prey biomass.
      ! This is specified by proportional_change=.true. The result of this is that the
      ! same specific prey loss rate (multiplied by distribution weights) will be applied over
      ! the predator depth range.
      call self%register_mapped_model_dependency(self%id_prey_int, 'prey', proportional_change=.true., domain=domain_surface)
   end subroutine

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_depth_integrated_predator), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: c, temp, prey_c, prey_n, prey_p, prey_s, w_int
      real(rk) :: ingestion_c, ingestion_n, ingestion_p, prey_loss_rate, p, net_growth
      integer  :: istate

      _SURFACE_LOOP_BEGIN_
         ! Get depth-integrated predator biomass
         _GET_SURFACE_(self%id_c, c)

         ! Depth-averaged environmental dependencies and prey concentrations
         _GET_SURFACE_(self%id_temp, temp)
         _GET_SURFACE_(self%id_prey_c, prey_c)
         _GET_SURFACE_(self%id_prey_n, prey_n)
         _GET_SURFACE_(self%id_prey_p, prey_p)

         ! Calculate ingested fluxes of different chemical elements
         ! Predator population growth will be based on the most limiting of these
         ingestion_c = self%clearance_rate * c * prey_c
         ingestion_n = self%clearance_rate * c * prey_n
         ingestion_p = self%clearance_rate * c * prey_p
         net_growth = min(ingestion_c, ingestion_n / NC, ingestion_p / PC) - self%mortality * c

         ! The specific loss rate of prey is the depth-integrated ingestion,
         ! divided by depth-integrated prey biomass, e.g., ingestion_c / prey_c_int.
         ! In turn, prey_c_int is related to depth-averaged prey as prey_c = prey_c_int / w_int,
         ! with w_int representing the depth-integral weights of the predator's vertical distibution.
         ! Thus, the specific loss rate is ingestion_c / (prey_c * w_int), which simplifies to
         ! clearance_rate * c / w_int (see expression for ingestion_c above)
         _GET_SURFACE_(self%id_w%integral, w_int)
         prey_loss_rate = self%clearance_rate * c / w_int

         ! Source term for predator
         _ADD_SURFACE_SOURCE_(self%id_c, net_growth)

         ! Apply the same specific loss rate of all state variables of the prey
         do istate = 1, size(self%id_prey_int%surface_state)
            _GET_SURFACE_(self%id_prey_int%surface_state(istate), p)
            _ADD_SURFACE_SOURCE_(self%id_prey_int%surface_state(istate), -prey_loss_rate * p)
         end do

         ! Send unused ingested matter and dead biomass to waste pools
         _ADD_SURFACE_SOURCE_(self%id_waste_c, ingestion_c - net_growth)
         _ADD_SURFACE_SOURCE_(self%id_waste_n, ingestion_n - net_growth * NC)
         _ADD_SURFACE_SOURCE_(self%id_waste_p, ingestion_p - net_growth * PC)

         ! Save diagnostics
         _SET_SURFACE_DIAGNOSTIC_(self%id_net_growth, net_growth * 86400.0_rk)
         _SET_SURFACE_DIAGNOSTIC_(self%id_prey_loss_rate, prey_loss_rate * 86400.0_rk)
         
      _SURFACE_LOOP_END_
   end subroutine

end module

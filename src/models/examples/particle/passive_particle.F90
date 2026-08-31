#include "fabm_driver.h"

module examples_particle_passive_particle

   ! A passive particle
   
   ! Two types are provided:
   ! * type_passive_particle_cs has constant C:N:P stoichoiometry.
   !   Its concentration is expressed in carbon units.
   !   Elemental ratios N:C and  P:C are configurable and default to Redfield.
   ! * type_passive_particle_vs has variable C:N:P stoichoiometry.

   use fabm_particle
   use fabm_types

   implicit none

   private

   type, extends(type_particle_model), public :: type_passive_particle_cs
      type (type_state_variable_id) :: id_c
   contains
      procedure :: initialize => initialize_cs
   end type

   type, extends(type_particle_model), public :: type_passive_particle_vs
      type (type_state_variable_id) :: id_c
      type (type_state_variable_id) :: id_n
      type (type_state_variable_id) :: id_p
   contains
      procedure :: initialize => initialize_vs
   end type

contains

   subroutine initialize_cs(self, configunit)
      class (type_passive_particle_cs), intent(inout), target :: self
      integer,                          intent(in)            :: configunit

      real(rk) :: NC, PC

      call self%register_state_variable(self%id_c, 'c', 'mmol C m-3', 'concentration')

      call self%get_parameter(NC, 'NC', 'mol N mol C-1', 'nitrogen : carbon ratio', default=16.0_rk/ 106.0_rk)
      call self%get_parameter(PC, 'PC', 'mol P mol C-1', 'phosphorus : carbon ratio', default=1.0_rk/ 106.0_rk)

      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c, scale_factor=NC)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=PC)
   end subroutine

   subroutine initialize_vs(self, configunit)
      class (type_passive_particle_vs), intent(inout), target :: self
      integer,                          intent(in)            :: configunit

      call self%register_state_variable(self%id_c, 'c', 'mmol C m-3', 'carbon')
      call self%register_state_variable(self%id_n, 'n', 'mmol N m-3', 'nitrogen')
      call self%register_state_variable(self%id_p, 'p', 'mmol P m-3', 'phosphorus')

      call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_c)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_n)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_p)
   end subroutine

end module

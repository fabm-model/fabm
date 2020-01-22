#include "fabm_driver.h"

module su_pools

  use fabm_types

  implicit none

  private

  type, extends(type_base_model), public :: type_su_pools
    ! Variables identifiers
    type(type_state_variable_id) :: id_N, id_P, id_C, id_S

    logical :: has_nitrogen, has_phosphorus, has_carbon, has_silicate
  contains
    ! Add model procedures here
    procedure :: initialize
  end type type_su_pools

contains

  subroutine initialize(self,configunit)                   
    class(type_su_pools), intent(inout), target :: self
    integer,              intent(in)            :: configunit

    ! Register variables in FABM
    call self%get_parameter(self%has_nitrogen,   'has_nitrogen',   '', 'whether this pool contains nitrogen',   default = .false.)
    if (self%has_nitrogen) call self%register_state_variable(self%id_N, 'N', 'ugN L-1', 'nitrogen')

    call self%get_parameter(self%has_phosphorus, 'has_phosphorus', '', 'whether this pool contains phosphorus', default = .false.)
    if (self%has_phosphorus) call self%register_state_variable(self%id_P, 'P', 'ugP L-1', 'phosphorus')

    call self%get_parameter(self%has_carbon,     'has_carbon',     '', 'whether this pool contains carbon',     default = .false.)
    if (self%has_carbon) call self%register_state_variable(self%id_C, 'C', 'ugC L-1', 'carbon')

    call self%get_parameter(self%has_silicate,   'has_silicate',   '', 'whether this pool contains silicate',   default = .false.)
    if (self%has_silicate) call self%register_state_variable(self%id_S, 'S', 'ugSi L-1', 'silicate')

    ! Register the contribution of variables to total N, P, C and S
    if (self%has_nitrogen) call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_N,  scale_factor=1._rk/14._rk)
    if (self%has_phosphorus) call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_P,  scale_factor=1._rk/31._rk)
    if (self%has_carbon) call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_C,  scale_factor=1._rk/12._rk)
    if (self%has_silicate) call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_S,  scale_factor=1._rk/60._rk)
  end subroutine initialize

end module su_pools

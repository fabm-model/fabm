#include "fabm_driver.h"

module su_npz

  use fabm_types

  implicit none

  private

  type,extends(type_base_model),public :: type_su_npz
    ! Variables identifiers
    type(type_state_variable_id) :: id_n, id_p, id_z
    type(type_diagnostic_variable_id) :: id_phy_u, id_zoo_u

    ! Parameters
    real (rk) :: AE_N
    real (rk) :: bas_res
    real (rk) :: phy_umax
    real (rk) :: ing_k
    real (rk) :: phy_k
    real (rk) :: SDA
    real (rk) :: thres
    real (rk) :: zoo_umax
  contains
    ! Add model procedures here
    procedure :: initialize
    procedure :: do
  end type

contains

  subroutine initialize(self,configunit)
    class(type_su_npz),intent(inout),target :: self
    integer,           intent(in)           :: configunit

    ! Register state variables
    call self%register_state_variable(self%id_n, 'n', 'ugN L-1', 'nitrogen')
    call self%register_state_variable(self%id_p, 'p', 'ugN L-1', 'phytoplankton')
    call self%register_state_variable(self%id_z, 'z', 'ugN L-1', 'zooplankton')

    ! Register the contribution of all state variables to total nitrogen
    call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n,scale_factor=1._rk/14._rk)
    call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_p,scale_factor=1._rk/14._rk)
    call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_z,scale_factor=1._rk/14._rk)

    ! Register diagnostic variables
    call self%register_diagnostic_variable(self%id_phy_u,'phy_u','gN (gN)-1 d-1','phytoplankton growth')
    call self%register_diagnostic_variable(self%id_zoo_u,'zoo_u','gN (gN)-1 d-1','zooplankton growth')

    ! Set parameter values
    call self%get_parameter(self%AE_N,'AE_N','dl','nitrogen assimilation efficiency',default=0.6_rk)
    call self%get_parameter(self%bas_res,'bas_res','gN (gN)-1 d-1','basal respiration rate',default=0.05_rk)
    call self%get_parameter(self%phy_umax,'phy_umax','gN (gN)-1 d-1','maximum phytoplankton growth rate',default=0.693_rk)
    call self%get_parameter(self%ing_k,'ing_k','ugN L-1','half-saturation constant for ingestion',default=28._rk)
    call self%get_parameter(self%phy_k,'phy_k','ugN L-1','half-saturation constant for nitrogen uptake',default=14._rk)
    call self%get_parameter(self%SDA,'SDA','dl','specific dynamic action',default=0.3_rk)
    call self%get_parameter(self%thres,'thres','ugN L-1','threshold value for predation',default=0.1_rk)
    call self%get_parameter(self%zoo_umax,'zoo_umax','gN (gN)-1 d-1','maximum zooplankton growth rate',default=1._rk)
  end subroutine initialize
  
  !add model subroutines here
  subroutine do(self, _ARGUMENTS_DO_)
    class (type_su_npz), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    !Local variables
    real(rk) :: n,p,z
    real(rk) :: d_n,d_p,d_z
    real(rk) :: phy_u, zoo_u
    real(rk) :: assi, ing, ing_max, phy_gro, reg, sys, zoo_ing, zoo_regen, zoo_void

    !If any spatial loops:
    _LOOP_BEGIN_

    !Retrieve current (local) state variable values
     _GET_(self%id_n,n)
     _GET_(self%id_p,p)
     _GET_(self%id_z,z)

    !Calculate auxiliaries
    ing_max   = (self%zoo_umax + self%bas_res)/(self%AE_N*(1._rk-self%SDA))
    if (p >= self%thres) then
       ing = ing_max*(p - self%thres)/(p - self%thres + self%ing_k)
    else
       ing = 0.0_rk
    end if
    assi      = ing*self%AE_N
    phy_u     = self%phy_umax*n/(n + self%phy_k)
    phy_gro   = phy_u*p
    reg       = self%bas_res + assi*self%SDA
    zoo_ing   = ing*z
    zoo_regen = reg*z
    zoo_u     = ing*self%AE_N*(1._rk - self%SDA) - self%bas_res
    zoo_void  = ing*(1._rk - self%AE_N)*z

    !Calculate temporal derivatives
    d_n = zoo_void + zoo_regen - phy_gro
    d_p = phy_gro - zoo_ing
    d_z = zoo_ing - zoo_void - zoo_regen

    !Provide temporal derivatives to FABM (must be in per second, hence division by 86400)
    _ADD_SOURCE_(self%id_n,d_n/86400)
    _ADD_SOURCE_(self%id_p,d_p/86400)
    _ADD_SOURCE_(self%id_z,d_z/86400)

    !Provide diagnostic variables to FABM
    _SET_DIAGNOSTIC_(self%id_phy_u,phy_u)
    _SET_DIAGNOSTIC_(self%id_zoo_u,zoo_u)
     
    !End spatial loops
    _LOOP_END_
  end subroutine do

end module su_npz


















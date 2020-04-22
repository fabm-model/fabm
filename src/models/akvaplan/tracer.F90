#include "fabm_driver.h"

module akvaplan_tracer

   ! A simple tracer model with support for sinking/floating and temperature-dependent decay.
   ! The temperature dependence of decay can be described with a Q10 formulation or an Arrhenius formulation.
   ! Copyright (C) 2016 - Akvaplan-niva

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_tracer
      ! The model's own variables
      type (type_state_variable_id) :: id_c ! tracer concentration

      ! Environmental dependencies
      type (type_dependency_id) :: id_T ! temperature

      ! Parameters
      real(rk) :: w
      real(rk) :: k
      integer :: temperature_dependence
      real(rk) :: T_ref
      real(rk) :: Q10
      real(rk) :: E_a
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type

   type (type_bulk_standard_variable),parameter :: non_water_volume_fraction = type_bulk_standard_variable(name='non_water_volume_fraction',units='-',      aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: non_water_density         = type_bulk_standard_variable(name='non_water_density',        units='kg m-3', aggregate_variable=.true.)

contains

   subroutine initialize(self,configunit)
      class (type_tracer),intent(inout),target :: self
      integer,            intent(in)           :: configunit

      real(rk)                        :: sp_vol, sp_dens
      logical                         :: conserved
      character(len=attribute_length) :: standard_name

      ! Obtain parameter values (convert all rate constants to s-1 to match FABM's internal time unit)
      call self%get_parameter(self%w, 'w', 'm d-1', 'sinking velocity',                default=0.0_rk, scale_factor=1.0_rk/86400)
      call self%get_parameter(self%k, 'k', 'd-1',   'decay rate (half-life = ln(2)/k)',default=0.0_rk, scale_factor=1.0_rk/86400)
      call self%get_parameter(self%temperature_dependence,'temperature_dependence', '', 'temperature dependence of decay (0: none, 1: Q10, 2: Arrhenius)', default=0, minimum=0, maximum=2)
      if (self%temperature_dependence/=0) call self%get_parameter(self%T_ref, 'T_ref', 'degree_C', 'temperature for which decay rate k is given')
      if (self%temperature_dependence==1) call self%get_parameter(self%Q10,   'Q10',   '-',        'Q10 temperature coefficient (1: no temperature dependence)', minimum=1.0_rk)
      if (self%temperature_dependence==2) call self%get_parameter(self%E_a,   'E_a',   'J mol-1',  'activation energy for decay (0: no temperature dependence)', minimum=0.0_rk)

      ! Register the model's own variables.
      call self%register_state_variable(self%id_c,'c','quantity m-3','concentration',initial_value=1.0_rk,vertical_movement=-self%w,minimum=0.0_rk)

      ! Register environmental dependencies.
      call self%register_dependency(self%id_T,standard_variables%temperature)

      ! Density hook based on specific volume and density of the tracer.
      call self%get_parameter(sp_vol,  'specific_volume', 'm3 quantity-1', 'specific volume', default=0.0_rk)
      call self%get_parameter(sp_dens, 'density','kg m-3', 'density (tracer mass/tracer volume)', default=0.0_rk)
      call self%add_to_aggregate_variable(non_water_volume_fraction,self%id_c,scale_factor=sp_vol)
      call self%add_to_aggregate_variable(non_water_density,        self%id_c,scale_factor=sp_vol*sp_dens) ! converting from kg tracer/tracer_volume to kg tracer/total_volume

      ! Optionally register the tracer as a "conserved quantity". This will prompt the hydrodynamic model to compute budgets across the entire domain.
      ! The name of the conserved quantity will be the short name of the model itself, with "_total" appended.
      call self%get_parameter(conserved,'conserved','','whether this is a conserved quantity (activates budget tracking in host)',default=.false.)
      if (conserved) then
         standard_name = get_safe_name(trim(self%get_path())//'_total')
         call self%add_to_aggregate_variable(type_bulk_standard_variable(name=standard_name(2:),units='quantity m-3',aggregate_variable=.true.,conserved=.true.),self%id_c)
      end if

   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_tracer),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk), parameter :: R = 8.3144598_rk    ! universal gas constant (J mol-1 K-1)
      real(rk), parameter :: Kelvin = 273.15_rk  ! offset between degrees Celsius and Kelvin

      real(rk) :: c, T, f

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Obtain model state and environment from FABM
         _GET_(self%id_c,c) ! tracer concentration
         _GET_(self%id_T,T) ! temperature (degree Celsius)

         ! Temperature-dependent scale factor (dimensionless) for decay rate
         if (self%temperature_dependence==1) then
            ! Q10 relationship (increase rate with a factor Q10 for every 10 degrees increase in temperature)
            f = self%Q10**((T-self%T_ref)/10)
         elseif (self%temperature_dependence==2) then
            ! Arrhenius relationship (temperature dependence of chemical reactions based on their activation energy E_a)
            f = exp(self%E_a/(R*(Kelvin+self%T_ref))-self%E_a/(R*(Kelvin+T)))
         else
            ! No temperature dependence
            f = 1.0_rk
         end if

         ! Send rate of change to FABM
         _ADD_SOURCE_(self%id_c,-self%k*f*c)

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

end module

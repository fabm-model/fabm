#include "fabm_driver.h"

module akvaplan_tracer_sed

   ! A simple tracer model with support for sinking/floating and temperature-dependent decay.
   ! The temperature dependence of decay can be described with a Q10 formulation or an Arrhenius formulation.
   ! Copyright (C) 2016 - Akvaplan-niva

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_tracer_sed
      ! The model's own variables
      type (type_state_variable_id)            :: id_c ! tracer concentration
      type (type_bottom_state_variable_id)     :: id_c_bot !

      ! Environmental dependencies
      type (type_dependency_id) :: id_T ! temperature
      type (type_bottom_dependency_id) :: id_shear !Bottom shear
      type (type_dependency_id) :: id_h ! layer thicknesis

      ! Parameters
      real(rk) :: w
      real(rk) :: k
      integer :: temperature_dependence
      logical :: do_sed
      integer :: resusp_meth
      real(rk) :: crt_shear,erate,bed_por
      real(rk) :: T_ref
      real(rk) :: Q10
      real(rk) :: E_a
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
   end type

   type (type_bulk_standard_variable),parameter :: non_water_volume_fraction = type_bulk_standard_variable(name='non_water_volume_fraction',units='-',      aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: non_water_density         = type_bulk_standard_variable(name='non_water_density',        units='kg m-3', aggregate_variable=.true.)
   
   !  private data members
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.0000000001_rk
   !real(rk),parameter :: max_dt=900._rk    !Maximum time step (s) to be expected - See Jorn's code   
   real(rk),parameter :: max_dt=0.64_rk     !OBS: Time step hardcoded from fvcom nml file (=EXTSTEP_SECONDS*ISPLIT)
contains

   subroutine initialize(self,configunit)
      class (type_tracer_sed),intent(inout),target :: self
      integer,            intent(in)           :: configunit

      real(rk)                        :: sp_vol, sp_dens
      logical                         :: conserved
      character(len=attribute_length) :: standard_name

      ! Obtain parameter values (convert all rate constants to s-1 to match FABM's internal time unit)
      call self%get_parameter(self%w, 'w', 'm d-1', 'sinking velocity',                default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%k, 'k', 'd-1',   'decay rate (half-life = ln(2)/k)',default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
      ! Do sedimentation, default is false
      call self%get_parameter(self%do_sed,'do_sed','','sedimentation switch (.true.: include sedimentation, .false.: disregard sedimentation)',default=.false.)
      ! Include resuspension, default = 0 (no resusp): build out gradually
      call self%get_parameter(self%resusp_meth,'resusp_meth','','resuspension method (0: no resuspension)', default=0, minimum=0,maximum=1)
      !critical shear stress (bottom) for resuspension
      call self%get_parameter(self%crt_shear,    'crt_shear',    'N m-2',             'critical shear stress',                                     default=0.005_rk)
      call self%get_parameter(self%erate,        'erate',        'kg m-2 s-1',        'bed erodibility constant',                                  default=0.0_rk)
      call self%get_parameter(self%bed_por,      'bed_por',      'm3 m-3',            'porosity [vol of voids/total vol]',                         default=0.0_rk)
     
      call self%get_parameter(self%temperature_dependence,'temperature_dependence', '', 'temperature dependence of decay (0: none, 1: Q10, 2: Arrhenius)', default=0, minimum=0, maximum=2)
      if (self%temperature_dependence/=0) call self%get_parameter(self%T_ref, 'T_ref', 'degree_C', 'temperature for which decay rate k is given')
      if (self%temperature_dependence==1) call self%get_parameter(self%Q10,   'Q10',   '-',        'Q10 temperature coefficient (1: no temperature dependence)', minimum=1.0_rk)
      if (self%temperature_dependence==2) call self%get_parameter(self%E_a,   'E_a',   'J mol-1',  'activation energy for decay (0: no temperature dependence)', minimum=0.0_rk)

      ! Register the model's own variables.
      call self%register_state_variable(self%id_c,'c','quantity m-3','concentration',initial_value=1.0_rk,vertical_movement=-self%w,minimum=0.0_rk)   
      call self%register_state_variable(self%id_c_bot,'c_bot','quantity m-2','concentration at bottom',minimum=0.0_rk)
  
      ! Register environmental dependencies.
      call self%register_dependency(self%id_T,standard_variables%temperature)
      call self%register_dependency(self%id_shear,  standard_variables%bottom_stress)
      call self%register_dependency(self%id_h,standard_variables%cell_thickness)

      ! Density hook based on specific volume and density of the tracer.
      call self%get_parameter(sp_vol,  'specific_volume', 'm3 quantity-1', 'specific volume', default=0.0_rk)
      call self%get_parameter(sp_dens, 'density','kg m-3', 'density (tracer mass/tracer volume)', default=0.0_rk)
      call self%add_to_aggregate_variable(non_water_volume_fraction,self%id_c,scale_factor=sp_vol)
      call self%add_to_aggregate_variable(non_water_volume_fraction,self%id_c_bot,scale_factor=sp_vol)
      call self%add_to_aggregate_variable(non_water_density,        self%id_c,scale_factor=sp_vol*sp_dens) ! converting from kg tracer/tracer_volume to kg tracer/total_volume
      call self%add_to_aggregate_variable(non_water_density,        self%id_c_bot,scale_factor=sp_vol*sp_dens)

      ! Optionally register the tracer as a "conserved quantity". This will prompt the hydrodynamic model to compute budgets across the entire domain.
      ! The name of the conserved quantity will be the short name of the model itself, with "_total" appended.
      call self%get_parameter(conserved,'conserved','','whether this is a conserved quantity (activates budget tracking in host)',default=.false.)
      if (conserved) then
         standard_name = get_safe_name(trim(self%get_path())//'_total')
         call self%add_to_aggregate_variable(type_bulk_standard_variable(name=standard_name(2:),units='quantity m-3',aggregate_variable=.true.,conserved=.true.),self%id_c)
      end if

   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_tracer_sed),intent(in) :: self
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
   
   !EOC
   
   !---------------------------------------------------------------------------------------------
   !BOP
   !
   !  !IROUTINE:
   !  
   !  !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   !
   !  !DESCRIPTION
   !
   !
   !  !INPUT PARAMETERS:
      class (type_tracer_sed),intent(in) :: self
         _DECLARE_ARGUMENTS_DO_BOTTOM_
   ! !REVISION HISTORY:
   !  Original author(s):
   !
   ! !LOCAL VARIABLES:
      real(rk) :: c
      real(rk) :: c_bot
      real(rk) :: shear,erosion,h
   !EOP
   !----------------------------------------------------------------------------------------------
   !BOC
      ! Enter spatial loops (if any)
         _HORIZONTAL_LOOP_BEGIN_
  
          ! Environment
          _GET_(self%id_c,c)
          _GET_(self%id_h,h)
          _GET_BOTTOM_(self%id_c_bot,c_bot)
          _GET_BOTTOM_(self%id_shear,shear)

          ! Sedimentation
          if (self%do_sed) then
            !_ADD_BOTTOM_FLUX_(self%id_c,-self%w*c)
            !_ADD_BOTTOM_SOURCE_(self%id_c_bot, +self%w*c)

                ! Resuspension/Erosion
                select case(self%resusp_meth)

                  case(0)
                  !Do nothing
                    erosion = 0.0_rk
                  case(1)
                  !resuspension as described in Ariathurai and Arulanandan (1978). 
                  !Also used in fvcom and ROMS (Warner et al. 2008)
                    
                    erosion = min(c_bot/max_dt, self%erate*(1.0_rk - self%bed_por)*max(shear/self%crt_shear - 1.0_rk, 0.0_rk))
                    !erosion = self%erate*(1.0_rk - self%bed_por)*max(shear/self%crt_shear - 1.0_rk, 0.0_rk)
                    
                    !write(*,*) 'self%erate ', self%erate
                    !write(*,*) 'self%bed_por', self%bed_por
                    !write(*,*) 'self%crt_shear', self%crt_shear 
                    !write(*,*) 'shear ', shear
                    !write(*,*) 'max(shear/self%crt_shear - 1.0_rk, 0.0_rk) ', max(shear/self%crt_shear - 1.0_rk, 0.0_rk) 
                    !write(*,*) 'erosion ', erosion
                end select

                
                
                !write(*,*) 'c_bot',c_bot
                !write(*,*) 'wtimesc',self%w*c
                !write(*,*) 'erosion',erosion
                !write(*,*) 'result',c_bot+self%w*c-erosion
                !write(*,*) 'id_c_bot',self%id_c_bot
                !write(*,*) 'h = ',h
                !write(*,*) 'max_dt = ',max_dt
                
                _ADD_BOTTOM_FLUX_(self%id_c,-min(h/max_dt,self%w)*c+erosion)
                _ADD_BOTTOM_SOURCE_(self%id_c_bot, +min(h/max_dt,self%w)*c-erosion)
                !_ADD_BOTTOM_FLUX_(self%id_c,-self%w*c+erosion)
                !_ADD_BOTTOM_SOURCE_(self%id_c_bot, +self%w*c-erosion)

          end if

         _HORIZONTAL_LOOP_END_
  
        end subroutine do_bottom
   !EOC
   
end module

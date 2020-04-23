#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_diatoms - diatoms lifestage model
!
! !INTERFACE:
   module fabm_uhh_diatoms
!
! !DESCRIPTION:
!
! The cold-water diatoms lifestage model
! from phD-thesis of Alexandra Warns (University of Hamburg)
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uhh_diatoms
!     Variable identifiers
      type (type_state_variable_id)         :: id_veg,id_res
      type (type_state_variable_id)         :: id_nitrate,id_ammonium
      type (type_state_variable_id)         :: id_phosphate,id_detritus,id_oxygen
      type (type_dependency_id)             :: id_par,id_temp
      type (type_global_dependency_id)      :: id_doy
      type (type_surface_state_variable_id) :: ssv_mean_growth
      type (type_diagnostic_variable_id)    :: id_gr, id_tc, id_llim, id_nlim, id_tlim, id_pard      

!     Model parameters
      real(rk) :: mumax
      real(rk) :: vmort
      real(rk) :: rmort
      real(rk) :: rkn
      real(rk) :: rdepo
      real(rk) :: alpha
      real(rk) :: w_dia
      real(rk) :: rkc
      real(rk) :: mdt
      real(rk) :: minimum_nitrate
      real(rk) :: trate_veg_res
      real(rk) :: trate_res_veg
      real(rk) :: sr
      real(rk) :: s2
      real(rk) :: s3
      logical  :: use_phosphate
      logical  :: use_ammonium
      logical  :: use_oxygen
      logical  :: use_mean_growth

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: check_surface_state
      procedure :: get_light_extinction
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------

   contains


   !> Initialise the diatoms model
   !
   !>  Here, the uhh_diatoms namelist is read and variables exported
   !!  by the model are registered with FABM.
   subroutine initialize(self,configunit)
   class (type_uhh_diatoms), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

   real(rk)     :: mumax= 0.4
   real(rk)     :: vmort=0.03
   real(rk)     :: rmort=0.03
   real(rk)     :: rkn=0.5
   real(rk)     :: rdepo=5.6e-6
   real(rk)     :: alpha=0.18
   real(rk)     :: w_dia=0.0
   real(rk)     :: rkc=0.07  
   real(rk)     :: sr=0.0625 
   real(rk)     :: s2=6.625
   real(rk)     :: s3=8.625 
   real(rk)     :: mdt=1800.
   real(rk)     :: minimum_nitrate=0.0_rk
   real(rk)     :: trate_veg_res     = 0.02
   real(rk)     :: trate_res_veg     = 0.02
   logical      :: use_mean_growth=.false.

   character(len=64)         :: phosphate_variable
   character(len=64)         :: ammonium_variable
   character(len=64)         :: nitrate_variable
   character(len=64)         :: detritus_variable
   character(len=64)         :: oxygen_variable
  namelist /uhh_diatoms/ mumax,vmort,rmort,rkn,rdepo,sr,         &
                      trate_veg_res,trate_res_veg,            &
                      w_dia,rkc,alpha,mdt,use_mean_growth,    &
                      ammonium_variable,nitrate_variable,     &
                      phosphate_variable, detritus_variable,  &
                      oxygen_variable, minimum_nitrate

   nitrate_variable = 'uhh_ergom_split_base_nit'
   ammonium_variable = 'uhh_ergom_split_base_amm'
   phosphate_variable = 'uhh_ergom_split_base_pho'
   detritus_variable = 'uhh_ergom_split_base_det'
   oxygen_variable = 'uhh_ergom_split_base_oxy'
   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_diatoms,err=99)

   ! set dependency switches
   self%use_phosphate = phosphate_variable /= ''
   self%use_ammonium  = ammonium_variable /= ''
   self%use_oxygen    = oxygen_variable /= ''

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%mumax, 'mumax', default=mumax,scale_factor=one_pr_day)
   call self%get_parameter(self%vmort, 'vmort', default=vmort,scale_factor=one_pr_day)
   call self%get_parameter(self%rmort, 'rmort', default=rmort,scale_factor=one_pr_day)
   call self%get_parameter(self%rkn,   'rkn',   default=rkn)
   call self%get_parameter(self%rdepo, 'rdepo', default=rdepo,scale_factor=one_pr_day)
   call self%get_parameter(self%alpha, 'alpha', default=alpha,scale_factor=one_pr_day)
   call self%get_parameter(self%w_dia, 'w_dia', default=w_dia,scale_factor=one_pr_day)
   call self%get_parameter(self%rkc,   'rkc',   default=rkc)
   call self%get_parameter(self%sr,    'sr',    default=sr)
   call self%get_parameter(self%s2,   's2',   default=s2)
   call self%get_parameter(self%s3,   's3',   default=s3) 
   call self%get_parameter(self%mdt,   'mdt', default=mdt)
   call self%get_parameter(self%minimum_nitrate, 'minimum_nitrate', default=minimum_nitrate)
   call self%get_parameter(self%use_mean_growth, 'use_mean_growth', default=use_mean_growth)
   call self%get_parameter(self%trate_veg_res, 'trate_veg_res', default=trate_veg_res,scale_factor=one_pr_day)
   call self%get_parameter(self%trate_res_veg, 'trate_res_veg', default=trate_res_veg,scale_factor=one_pr_day)

   ! Register state variables
   call self%register_state_variable(self%id_veg,'veg', &
         'mmol n/m**3','diatoms vegetative biomass',  &
         minimum=0.0e-7_rk,vertical_movement=self%w_dia)

   call self%register_state_variable(self%id_res,'res', &
         'mmol n/m**3','vegetatives resting biomass', &
         minimum=0.0e-7_rk,vertical_movement=self%w_dia)

   if (self%use_mean_growth) &
         call self%register_surface_state_variable(self%ssv_mean_growth,'meangrowth', &
         '1/s','mean growth rate', initial_value=0.0_rk)

   ! Register dependencies on external standard variables
   if (self%use_ammonium) &
     call self%register_state_dependency(self%id_ammonium, 'ammonium_target', 'mmol/m**3','ammonium source')
   call self%register_state_dependency(self%id_nitrate, 'nitrate_target', 'mmol/m**3','nitrate source')
   if (self%use_phosphate) &
     call self%register_state_dependency(self%id_phosphate, 'phosphate_target',  'mmol/m**3','phosphate source')

   
   ! Register external state dependencies
   call self%register_state_dependency(self%id_detritus, 'mortality_target','mmol/m**3','sink for dead matter')
   if (self%use_oxygen) &
     call self%register_state_dependency(self%id_oxygen,   'oxygen_target'   ,'mmol-O2/m**3','dissolved oxygen pool')

   if (self%use_ammonium) &
     call self%request_coupling(self%id_ammonium, ammonium_variable)
   call self%request_coupling(self%id_nitrate, nitrate_variable)
   if (self%use_phosphate) &
     call self%request_coupling(self%id_phosphate, phosphate_variable)
   call self%request_coupling(self%id_detritus,detritus_variable)
   if (self%use_oxygen) &
     call self%request_coupling(self%id_oxygen, oxygen_variable)
   
   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_gr,'gr_veg','1/d', &
      'relative growth rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tc,'trans_ctrl','', &
      'transition from vegetative diatoms to spores', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_llim,'light_limitation','', &
      'light limitation', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_nlim,'nutrient_limitation','', &
      'nutrient limitation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tlim,'temp_limitation','', &
      'temperature limitation', output=output_instantaneous)      
   call self%register_diagnostic_variable(self%id_pard,'par','W/m**2', &
      'photosynthetically active radiation', output=output_instantaneous)    
      
   return

99 call self%fatal_error('fabm_uhh_diatoms','Error reading namelist uhh_diatoms')

   end subroutine initialize



   !> Right hand sides of diatoms model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uhh_diatoms), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

   real(rk) :: ni,am,d,po,par,temp,nut
   real(rk) :: veg,res
   real(rk) :: tlim,llim,nlim,gr_veg
   real(rk) :: mean_growth
   real(rk) :: tau_res_veg,tau_veg_res
   real(rk) :: doy
   real(rk) :: trans_ctrl

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_veg,veg)       ! vegetatives biomass
   _GET_(self%id_res,res)       ! resting cells biomass
   _GET_(self%id_nitrate,ni)    ! nitrate
   _GET_(self%id_detritus,d)    ! detritus
   if (self%use_ammonium) then
     _GET_(self%id_ammonium,am) ! ammonium
   else
     am=0.0_rk
   end if
   nut = ni + am

   if (self%use_phosphate) then
     _GET_(self%id_phosphate,po) ! phosphate
   else
     po = nut/16.0_rk
   end if

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)              ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)            ! local temperature
   _GET_GLOBAL_(self%id_doy,doy)       ! day of year

   ! vegetatives growth
   llim=1.0_rk - exp( -(self%alpha * par) / self%mumax)
   nlim=nut/(nut+self%rkn)
   tlim=0.5_rk *(tanh(0.5_rk * (temp - 2.5_rk)) - (tanh(0.65_rk *(temp - 15.0_rk)) ))
   gr_veg= self%mumax*llim*nlim*tlim

   ! either use instantaneous values or averaged growth rate
   if (self%use_mean_growth) then
     _GET_HORIZONTAL_(self%ssv_mean_growth,mean_growth)
   else
     mean_growth = gr_veg
   end if

   ! lifestage fluxes
   tau_res_veg=0.0_rk
   tau_veg_res=0.0_rk
   trans_ctrl = 0.5_rk * (1.0_rk - tanh(200.0_rk * secs_pr_day * &
                (mean_growth - (0.15_rk * one_pr_day))))
   ! veg -> res
   tau_veg_res=self%trate_veg_res * trans_ctrl
   ! res -> veg
   tau_res_veg=self%trate_res_veg * (1.0_rk - trans_ctrl)

   ! Set temporal derivatives
   _ADD_SOURCE_(self%id_veg,veg*(gr_veg - self%vmort) - veg*tau_veg_res + res*tau_res_veg)
   _ADD_SOURCE_(self%id_res,veg*tau_veg_res - res*(self%rmort + tau_res_veg))

   ni = max(ni, self%minimum_nitrate)
   
   ! external nutrients
   _ADD_SOURCE_(self%id_nitrate,-veg*gr_veg * ni/(ni+am))
   if (self%use_ammonium) then
     _ADD_SOURCE_(self%id_ammonium,-veg*gr_veg * am/(ni+am))
   end if
   _ADD_SOURCE_(self%id_detritus,veg*self%vmort + res*self%rmort)
   if (self%use_phosphate) then
     _ADD_SOURCE_(self%id_phosphate, -self%sr *veg*gr_veg)
   end if
   ! add oxygen dynamics
   _ADD_SOURCE_(self%id_oxygen, (self%s2 *am/(am+ni) + self%s3* ni/(am+ni))*veg*gr_veg)
   
   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_gr,gr_veg*secs_pr_day)   
   _SET_DIAGNOSTIC_(self%id_tc,trans_ctrl)
   _SET_DIAGNOSTIC_(self%id_llim,llim)
   _SET_DIAGNOSTIC_(self%id_nlim,nlim)
   _SET_DIAGNOSTIC_(self%id_tlim,tlim)
   _SET_DIAGNOSTIC_(self%id_pard,par)
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

   
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_uhh_diatoms), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
   real(rk),save :: res

   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_res,res)       ! resting cells biomass
   
   _ADD_BOTTOM_FLUX_(self%id_res,-self%rdepo*res*res)
   

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
   

   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
   class (type_uhh_diatoms), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
   
   real(rk)                     :: veg,res

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_veg,veg)       ! vegetatives biomass
   _GET_(self%id_res,res)       ! resting cells biomass

   ! Self-shading
   _SET_EXTINCTION_(self%rkc*(veg+res))

   _LOOP_END_

   end subroutine get_light_extinction


   subroutine check_surface_state(self,_ARGUMENTS_CHECK_SURFACE_STATE_)
   class (type_uhh_diatoms), intent(in) :: self
   _DECLARE_ARGUMENTS_CHECK_SURFACE_STATE_
   real(rk)            :: par,temp,nut,mean_growth,gr_veg,ni,am
   real(rk)            :: llim,nlim,tlim
   ! set averaging window (0.5 * original averaging window, because
   ! of dampening, recursive filter)
   real(rk), parameter :: window = 2.5_rk*86400.0_rk
   real(rk), parameter :: one_pr_window = 1.0_rk/window

   if (.not.self%use_mean_growth) return

   _HORIZONTAL_LOOP_BEGIN_
   _GET_HORIZONTAL_(self%ssv_mean_growth, mean_growth) ! mean growth

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)       ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)     ! local temperature
   _GET_(self%id_nitrate,ni)    ! nitrate
   if (self%use_ammonium) then
     _GET_(self%id_ammonium,am) ! ammonium
   else
     am=0.0_rk
   end if
   nut = ni + am

   ! vegetatives growth
   llim=1.0_rk - exp( -(self%alpha * par) / self%mumax)
   nlim=nut/(nut+self%rkn)
   tlim=0.5_rk *(tanh(0.5_rk * (temp - 2.5_rk)) - (tanh(0.65_rk *(temp - 15.0_rk)) ))
   gr_veg= self%mumax*llim*nlim*tlim
 
   mean_growth = one_pr_window*((window - self%mdt)*mean_growth + self%mdt*gr_veg)
   _SET_HORIZONTAL_(self%ssv_mean_growth,mean_growth)
   _HORIZONTAL_LOOP_END_

   end subroutine check_surface_state


   end module fabm_uhh_diatoms


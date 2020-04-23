#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_clcbase - CLC cyanobacteria lifestage base model
!
! !INTERFACE:
   module fabm_uhh_clcbase
!
! !DESCRIPTION:
!
! base model of a cyanobacteria lifestage within a cyanobacteria lifecycle.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uhh_clcbase
!     Variable identifiers
      type (type_state_variable_id)        :: id_c,id_s,id_g
      type (type_state_variable_id)        :: id_nitrate,id_ammonium
      type (type_state_variable_id)        :: id_phosphate,id_detritus,id_oxygen
      type (type_dependency_id)            :: id_par,id_temp, id_pml_carbonate_pH
      type (type_dependency_id)            :: id_salt
      type (type_state_variable_id)        :: id_nextc,id_nexts,id_nextg
      type (type_horizontal_dependency_id) :: id_I_0
      type (type_diagnostic_variable_id)   :: id_diagcflux,id_diaggflux,id_diagsflux
      type (type_diagnostic_variable_id)   :: id_diage,id_diagq
      type (type_diagnostic_variable_id)   :: id_up,id_gr,id_lc,id_nf,id_f_pH
      type (type_diagnostic_variable_id)   :: id_nf_mean
      type (type_diagnostic_variable_id)   :: id_primprod
      type (type_diagnostic_variable_id)   :: id_netprimprod 
      type (type_diagnostic_variable_id)   :: id_nfix 
      type (type_horizontal_diagnostic_variable_id) :: id_sediment_depo

!     Model parameters
      real(rk) :: sr
      real(rk) :: cya0
      real(rk) :: alpha
      real(rk) :: kc
      real(rk) :: s2
      real(rk) :: s3
      real(rk) :: mort
      real(rk) :: w
      real(rk) :: Qc
      real(rk) :: E_sc  ! energy storage capacity Emax in Hense&Beckmann 2006
      real(rk) :: kN
      real(rk) :: tscale
      real(rk) :: trange
      real(rk) :: tfc_c=25.0_rk
      real(rk) :: theta_e
      real(rk) :: theta_q
      real(rk) :: rmatscale
      real(rk) :: omega0
      real(rk) :: scale=8.0_rk
      real(rk) :: fcy2
      logical  :: n_fixation, lifecycling
      logical  :: use_ph
      real(rk) :: fpH_const
      real(rk) :: m
      real(rk) :: uptake_factor
      real(rk) :: growth_factor
      real(rk) :: lightcapture_factor
      real(rk) :: e_max,e_min,q_max,q_min
      real(rk) :: Sflux_per_Cflux, Gflux_per_Cflux
      real(rk) :: depo
      real(rk) :: minimum_C
      real(rk) :: minimum_nitrate
      character(len=64) :: next_in_cycle=''
      character(len=64) :: lifestage_name=''
      logical  :: use_phosphate
      logical  :: use_ammonium
      logical  :: use_oxygen
      logical  :: use_salinity_dependence

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
      procedure :: do_bottom
      procedure :: do_bottom_ppdd
      procedure :: get_light_extinction
      procedure :: init_lifestage => base_init_lifestage
      procedure :: calculate_lifecycle_flux => base_calculate_lifecycle_flux
   end type
!EOP
!-----------------------------------------------------------------------

   contains


   subroutine base_init_lifestage(self)
     class(type_uhh_clcbase) :: self
   end subroutine

   subroutine base_calculate_lifecycle_flux(self,e,q,c_flux,s_flux,g_flux)
     class(type_uhh_clcbase) :: self
     real(rk) :: e,q,c_flux,s_flux,g_flux
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the CLC model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the uhh_clc namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_uhh_clcbase), intent(inout), target :: self
   integer,                        intent(in)          :: configunit
!
! !LOCAL VARIABLES:
   real(rk)                  :: background_concentration
   real(rk)                  :: w
   real(rk)                  :: kc
   real(rk)                  :: mortality_rate
   real(rk)                  :: alpha
   real(rk)                  :: sr,s2,s3
   real(rk)                  :: Qc
   real(rk)                  :: kN
   real(rk)                  :: e_max,e_min,q_max,q_min
   real(rk)                  :: Sflux_per_Cflux, Gflux_per_Cflux
   real(rk)                  :: e_storage_capacity
   real(rk)                  :: tscale
   real(rk)                  :: trange
   real(rk)                  :: theta_e
   real(rk)                  :: theta_q
   real(rk)                  :: omega0
   real(rk)                  :: rmatscale
   real(rk)                  :: scale=8.0_rk
   logical                   :: n_fixation
   logical                   :: use_ph
   logical                   :: use_salinity_dependence
   real(rk)                  :: fpH_const=1.0_rk
   real(rk)                  :: m
   real(rk)                  :: fcy2
   real(rk)                  :: depo
   real(rk)                  :: uptake_factor
   real(rk)                  :: growth_factor
   real(rk)                  :: lightcapture_factor
   character(len=64)         :: next_in_cycle
   character(len=64)         :: phosphate_variable
   character(len=64)         :: ammonium_variable
   character(len=64)         :: nitrate_variable
   character(len=64)         :: detritus_variable
   character(len=64)         :: oxygen_variable
   character(len=64)         :: lifestage_name
   real(rk)                  :: minimum_C=0.0_rk
   real(rk)                  :: minimum_nitrate=0.0_rk

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   namelist /uhh_clc/ background_concentration, &
             w, kc, sr, s2, s3, rmatscale, omega0, alpha, fcy2, &
             Qc,kN,tscale,trange,theta_e,theta_q,scale, e_storage_capacity, &
             e_max,e_min,q_max,q_min,Sflux_per_Cflux, Gflux_per_Cflux, &
             mortality_rate, next_in_cycle, n_fixation, m, depo, &
             uptake_factor, growth_factor, lightcapture_factor, &
             ammonium_variable,nitrate_variable, use_ph, &
             use_salinity_dependence,                    &
             fpH_const, phosphate_variable, detritus_variable, &
             oxygen_variable, lifestage_name, minimum_C, minimum_nitrate
!EOP
!-----------------------------------------------------------------------
!BOC

   background_concentration = 0.0225_rk
   w         = 0.1_rk !sinking rate
   kc        = 0.03_rk !Self-shading kd
   sr        = 0.0625_rk
   s2        = 6.625_rk
   s3        = 8.625_rk 
   Qc        = 5.5_rk
   kN        = 0.3_rk
   e_max     = 1.0_rk !0~1
   e_min     = 0.0_rk !0~1
   q_max     = 1.0_rk !0~1
   q_min     = 0.0_rk !0~1
   Sflux_per_Cflux = 1.0_rk !0~1
   Gflux_per_Cflux = 1.0_rk !0~1
   e_storage_capacity = 5.5_rk
   tscale    = 2.5e-6_rk
   trange    = 0.05_rk
   theta_e   = 0.0_rk
   theta_q   = 0.0_rk
   scale     = 8.0_rk !omega
   rmatscale = 45.0_rk
   omega0    = 3.7e-6_rk
   n_fixation=.false.
   m         = 3.0_rk !grhs
   fcy2      = 4.0_rk !growth
   depo      = 0.0_rk !sedment deposition rate
   uptake_factor = 0.0_rk
   growth_factor = 0.0_rk
   lightcapture_factor = 0.0_rk
   mortality_rate = 0.02_rk ! mort [1/d]
   use_salinity_dependence=.false.

   nitrate_variable = 'uhh_ergom_split_base_nit'
   ammonium_variable = 'uhh_ergom_split_base_amm'
   phosphate_variable = 'uhh_ergom_split_base_pho'
   detritus_variable = 'uhh_ergom_split_base_det'
   oxygen_variable = 'uhh_ergom_split_base_oxy'
   lifestage_name = ''
   next_in_cycle = ''

   !use_ph=.true.
   use_ph=.false.
   fpH_const=1.0_rk

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_clc,err=99)

   ! set dependency switches
   self%use_phosphate = phosphate_variable /= ''
   self%use_ammonium  = ammonium_variable /= ''
   self%use_oxygen    = oxygen_variable /= ''
   self%use_ph        = use_ph

   !! call lifestage specific initialisation
   call self%init_lifestage()

   if (next_in_cycle /= '') self%next_in_cycle=next_in_cycle
   if (lifestage_name /= '') self%lifestage_name=lifestage_name
   self%lifecycling = self%next_in_cycle /= ''

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%cya0, 'background_concentration', default=background_concentration)
   call self%get_parameter(self%kc,     'kc',     default=kc)
   call self%get_parameter(self%alpha,  'alpha',  default=alpha)
   call self%get_parameter(self%sr,     'sr',     default=sr)
   call self%get_parameter(self%s2,     's2',     default=s2)
   call self%get_parameter(self%s3,     's3',     default=s3)
   call self%get_parameter(self%e_max,  'e_max',  default=e_max)
   call self%get_parameter(self%e_min,  'e_min',  default=e_min)
   call self%get_parameter(self%q_max,  'q_max',  default=q_max)
   call self%get_parameter(self%q_min,  'q_min',  default=q_min)
   call self%get_parameter(self%scale,  'scale',  default=scale)
   call self%get_parameter(self%use_salinity_dependence,  'use_salinity_dependence',  default=use_salinity_dependence)

   ! set default S,G flux scales depending on E,Q range:
   call self%get_parameter(self%Sflux_per_Cflux, 'Sflux_per_Cflux', default=Sflux_per_Cflux)
   call self%get_parameter(self%Gflux_per_Cflux, 'Gflux_per_Cflux', default=Gflux_per_Cflux)

   call self%get_parameter(self%e_sc, 'e_storage_capacity', default=e_storage_capacity)
   call self%get_parameter(self%tscale, 'tscale', default=tscale)
   call self%get_parameter(self%trange, 'trange', default=trange)
   call self%get_parameter(self%theta_e,'theta_e',default=theta_e)
   call self%get_parameter(self%theta_q,'theta_q',default=theta_q)
   call self%get_parameter(self%kN,     'kN',     default=kN)
   call self%get_parameter(self%Qc,     'Qc',     default=Qc)
   call self%get_parameter(self%m,      'm',      default=m)
   call self%get_parameter(self%fcy2,   'fcy2',   default=fcy2)
   call self%get_parameter(self%depo,   'depo',   default=depo)
   call self%get_parameter(self%n_fixation, 'n_fixation', default=n_fixation)
   call self%get_parameter(self%uptake_factor, 'uptake_factor', default=uptake_factor)
   call self%get_parameter(self%growth_factor, 'growth_factor', default=growth_factor)
   call self%get_parameter(self%lightcapture_factor, 'lightcapture_factor', default=lightcapture_factor)
   call self%get_parameter(self%omega0, 'omega0', default=omega0)
   call self%get_parameter(self%rmatscale, 'rmatscale', default=rmatscale)
   call self%get_parameter(self%w,      'w',      default=w,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%mort,   'mortality_rate', default=mortality_rate,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fpH_const,'fpH_const', default=fpH_const)
   call self%get_parameter(self%minimum_C,'minimum_C', default=minimum_C)
   call self%get_parameter(self%minimum_nitrate,'minimum_nitrate', default=minimum_nitrate)
   call self%get_parameter(self%lifestage_name,'lifestage_name', default=lifestage_name)
   call self%get_parameter(self%next_in_cycle,'next_in_cycle', default=next_in_cycle)
   call self%get_parameter(self%use_ph,'use_ph', default=use_ph)   
   ! Register state variables

   call self%register_state_variable(self%id_c,'C', &
         'mmol n/m**3',trim(self%lifestage_name)//' biomass',  &
         minimum=minimum_C,vertical_movement=self%w)

   call self%register_state_variable(self%id_g,'G', &
         'mmol n/m**3',trim(self%lifestage_name)//' energy', &
         minimum=0.3_rk*minimum_C,vertical_movement=self%w)

    call self%register_state_variable(self%id_s,'S', &
         'mmol n/m**3',trim(self%lifestage_name)//' quota',  &
         minimum=0.3_rk*minimum_C,vertical_movement=self%w)
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

   ! Register environmental dependencies
   if (self%lifecycling) then
      call self%register_state_dependency(self%id_nextc,'next_C','mmol-N/m**3','next clc biomass')
      call self%register_state_dependency(self%id_nextg,'next_G','mmol-N/m**3','next clc energy')
      call self%register_state_dependency(self%id_nexts,'next_S','mmol-N/m**3','next clc quota')
      call self%request_coupling(self%id_nextc,'uhh_clc'//self%next_in_cycle(1:3)//'_C')
      call self%request_coupling(self%id_nextg,'uhh_clc'//self%next_in_cycle(1:3)//'_G')
      call self%request_coupling(self%id_nexts,'uhh_clc'//self%next_in_cycle(1:3)//'_S')
   end if
 
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
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
   if (self%use_ph) then
     call self%register_dependency(self%id_pml_carbonate_pH,standard_variables%ph_reported_on_total_scale)
   end if
   call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_diagcflux,'Cflux','mmol-N/m**3/d', &
         'C-flux into next lifestage', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_diaggflux,'Gflux','mmol-N/m**3/d', &
         'G-flux into next lifestage', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_diagsflux,'Sflux','mmol-N/m**3/d', &
         'S-flux into next lifestage', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_diage,'E','1/1', &
         'energy status', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_diagq,'Q','1/1', &
         'internal quota', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_gr,'growth','1/d', &
         'relative growth rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_up,'uptake','1/d', &
         'relative uptake rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_lc,'lightcapture','1/d', &
         'relative light capturing rate', output=output_instantaneous)
   if (self%n_fixation) then
     call self%register_diagnostic_variable(self%id_nf_mean,'nfixation_mean','1/d', &
         'relative time-integrated fixation rate', output=output_time_step_averaged)
     call self%register_diagnostic_variable(self%id_nf,'nfixation','1/d', &
         'relative fixation rate', output=output_instantaneous)
   end if
   if (self%use_ph) call self%register_diagnostic_variable(self%id_f_pH,'f_pH','1/1', &
         'pH dependency factor', output=output_instantaneous)

   call self%register_diagnostic_variable(self%id_primprod,'primprod','mmol-N/m**3/d', &
         'primary production rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_netprimprod,'netprimprod','mmol-N/m**3/d', &
         'net primary production rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_nfix,'nfix','mmol-N/m**3/d', &
         'N-fixation rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_sediment_depo,'seddepo','mmol-N/m**2/d', &
         'sed-deposition rate', output=output_instantaneous, source=source_do_bottom)      
   return

99 call self%fatal_error('fabm_uhh_clc','Error reading namelist uhh_clc')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_uhh_clcbase), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: ni,am,c,g,s,d,po,par,I_0,temp,e,q,ep,fqc,nut
   real(rk)                   :: uptake,light_capture,s_flux,g_flux,c_flux
   real(rk)                   :: growth,fixation,grhs
   real(rk)                   :: salt,sigma_salt
   real(rk)                   :: pHval, f_ph
   real(rk)                   :: q_min,Sflux_per_Cflux          
   real(rk)                   :: sigma_e,sigma_n,sigma_l,sigma_q,omega
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
   real(rk), parameter        :: c0 = 1.0e-8_rk
   real(rk), parameter        :: g0 = 5.0e-10_rk
   real(rk), parameter        :: s0 = 5.0e-10_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_c,c)                ! biomass
   _GET_(self%id_g,g)                ! energy
   _GET_(self%id_s,s)                ! quota
   _GET_(self%id_nitrate,ni)         ! nitrate
   _GET_(self%id_detritus,d)         ! detritus
   
   if (self%use_ammonium) then
     _GET_(self%id_ammonium,am)      ! ammonium
   else
     am=0.0_rk
   end if
   nut = ni + am

   if (self%use_phosphate) then
     _GET_(self%id_phosphate,po)     ! phosphate
   else
     po = nut/16.0_rk
   end if

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)           ! local temperature
   _GET_(self%id_salt,salt)           ! local salinity
   if (self%use_ph) then
     _GET_(self%id_pml_carbonate_pH,pHval)           ! local pH
     ! pH dependency
     f_pH = exp(-((pHval - 8.1_rk)**2)/(0.6_rk**2))
   else
     f_pH = self%fpH_const
   end if
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation
   

   !RH: nan for t==12.0_rk !
   omega = self%omega0 * (1.0_rk/self%rmatscale + &
     1.0_rk/(0.25_rk+exp(3.0_rk/(temp-12.0_rk)-0.5_rk)+exp(-500.0_rk/(temp-12.0_rk)+self%tfc_c)))
   ! no growth for salinities below 5 and above 18
   if (self%use_salinity_dependence) then
     sigma_salt = max(min(1.0_rk,6.5_rk-abs(salt-11.5_rk)),0.0_rk)
   else
     sigma_salt=1.0_rk
   end if

   sigma_l = omega * self%alpha * par / sqrt(omega**2 + self%alpha**2 * par**2)
   sigma_n = nut / (nut + self%kN)
   ! apply salinity limitation to omega to be active for 
   ! growth, uptake and N-fixation
   omega = omega * sigma_salt

   e = 0.0_rk
   q = 0.0_rk
   if (c > 0.0_rk) then
   !RH: maybe relax to e0,q0 for low g,s,c (e=(g+g0)/(c+c0)
     e = g/(c+c0)
     q = s/(c+c0)
     !e = (g+g0)/(c+c0)
     !q = (s+s0)/(c+c0)
   endif
   ! previous "emax" is replaced by self%e_sc (in namelist: E_storage_capacity)
   ! previous "er" is replaced by sigma_e (see Hense&Beckmann(2006))
   ! previous "fq" is replaced by sigma_q (see Hense&Beckmann(2006))
   ep       = max(0.0_rk,tanh(self%scale-self%scale*e/self%e_sc))
   sigma_e  = tanh(           self%scale*e/self%e_sc)
   fqc      = max(0.0_rk,tanh(self%scale-self%scale*q/self%Qc))
   sigma_q  = tanh(           self%scale*q/self%Qc)

   ! specific rates 
   light_capture = f_pH * self%lightcapture_factor * sigma_l * ep
   uptake = f_pH * self%uptake_factor * omega * sigma_n * sigma_e * fqc
   growth = f_pH * self%growth_factor * omega * sigma_e * sigma_q
   ! fixation_factor is the energy consumption factor for fixation = 3 for heterocysts
   fixation = 0.0_rk
   if (self%n_fixation) then
     fixation = f_pH * omega * sigma_e
     growth = self%fcy2 * growth
   end if

   c_flux=0.0_rk
   g_flux=0.0_rk
   s_flux=0.0_rk

   if (self%lifecycling) call self%calculate_lifecycle_flux(e,q,c_flux,s_flux,g_flux)

     ! Set temporal derivatives in next lifecycle
   _ADD_SOURCE_(self%id_nextc, c*c_flux)
   _ADD_SOURCE_(self%id_nexts, s*s_flux)
   _ADD_SOURCE_(self%id_nextg, g*g_flux)

   ! Set temporal derivatives
   _ADD_SOURCE_(self%id_c,c*(growth + self%growth_factor*fixation - self%mort - c_flux))
   _ADD_SOURCE_(self%id_s,c*(uptake - growth) - s*(self%mort + s_flux))
    grhs = c*(light_capture - uptake - &
           (self%m + self%growth_factor)*fixation - growth) - &
           g*(self%mort + g_flux)
   _ADD_SOURCE_(self%id_g,grhs)

   ni = max(ni, self%minimum_nitrate)
   
   ! external nutrients
   _ADD_SOURCE_(self%id_nitrate,-c*uptake * ni/(ni+am))
   if (self%use_ammonium) then
     _ADD_SOURCE_(self%id_ammonium,-c*uptake * am/(ni+am))
   end if
   _ADD_SOURCE_(self%id_detritus,c*self%mort+s*self%mort)
   if (self%use_phosphate) then
      _ADD_SOURCE_(self%id_phosphate, -self%sr*c*(uptake + self%growth_factor*fixation))
   end if
   ! add oxygen dynamics
   _ADD_SOURCE_(self%id_oxygen, self%s3*c* ni/(am+ni)*(growth + self%growth_factor*fixation))  

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_diagcflux,c*c_flux*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_diaggflux,g*g_flux*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_diagsflux,s*s_flux*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_diage,e)
   _SET_DIAGNOSTIC_(self%id_diagq,q)
   _SET_DIAGNOSTIC_(self%id_gr,growth*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_up,uptake*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_lc,light_capture*86400.0_rk)
   if (self%n_fixation) then
     _SET_DIAGNOSTIC_(self%id_nf,self%growth_factor*fixation*86400.0_rk)
     _SET_DIAGNOSTIC_(self%id_nf_mean,self%growth_factor*fixation*86400.0_rk)
   end if
   if (self%use_ph) then
     _SET_DIAGNOSTIC_(self%id_f_pH,f_pH)
   end if
   _SET_DIAGNOSTIC_(self%id_primprod,c*(growth + self%growth_factor*fixation)*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_netprimprod,c*(growth + self%growth_factor*fixation - self%mort - c_flux)*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_nfix,c*(self%growth_factor*fixation)*secs_pr_day)  

! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_uhh_clcbase), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                     :: c
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_c,c) ! biomass concentration

   ! Self-shading
   _SET_EXTINCTION_(self%kc*c)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

   subroutine do_ppdd(self,_ARGUMENTS_DO_PPDD_)
   class (type_uhh_clcbase),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_

   real(rk)                   :: ni,am,c,g,s,po,par,I_0,temp,e,q,ep,fqc,nut
   real(rk)                   :: uptake,light_capture,s_flux,g_flux,c_flux
   real(rk)                   :: growth,fixation
   real(rk)                   :: pHval, f_ph
   real(rk)                   :: salt, sigma_salt
   real(rk)                   :: sigma_e,sigma_n,sigma_l,sigma_q,omega
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
   real(rk), parameter        :: c0 = 1.0e-8_rk
   real(rk), parameter        :: g0 = 5.0e-10_rk
   real(rk), parameter        :: s0 = 5.0e-10_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_c,c)                ! biomass
   _GET_(self%id_g,g)                ! energy
   _GET_(self%id_s,s)                ! quota
   _GET_(self%id_nitrate,ni)         ! nitrate
   if (self%use_ammonium) then
     _GET_(self%id_ammonium,am)      ! ammonium
   else
     am=0.0_rk
   end if
   nut = ni + am

   if (self%use_phosphate) then
     _GET_(self%id_phosphate,po)     ! phosphate
   else
     po = nut/16.0_rk
   end if

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)           ! local temperature
   _GET_(self%id_salt,salt)           ! local salinity
   if (self%use_ph) then
     _GET_(self%id_pml_carbonate_pH,pHval)           ! local pH
     ! pH dependency
     f_pH = exp(-((pHval - 8.1_rk)**2)/(0.6_rk**2))
   else
     f_pH = self%fpH_const
   end if
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation


   !RH: nan for t==12.0_rk !
   omega = f_pH * self%omega0 * (1.0_rk/self%rmatscale + &
     1.0_rk/(0.25_rk+exp(3.0_rk/(temp-12.0_rk)-0.5_rk)+exp(-500.0_rk/(temp-12.0_rk)+self%tfc_c)))
   ! no growth for salinities below 5 and above 18
   if (self%use_salinity_dependence) then
     sigma_salt = max(min(1.0_rk,6.5_rk-abs(salt-11.5_rk)),0.0_rk)
   else
     sigma_salt=1.0_rk
   end if

   sigma_l = omega * self%alpha * par / sqrt(omega**2 + self%alpha**2 * par**2)
   sigma_n = nut / (nut + self%kN)

   e = 0.0_rk
   q = 0.0_rk
   if (c > 0.0_rk) then
   !RH: maybe relax to e0,q0 for low g,s,c (e=(g+g0)/(c+c0)
     e = g/c
     q = s/c
     !e = g/(c+c0)
     !q = s/(c+c0)
     !e = (g+g0)/(c+c0)
     !q = (s+s0)/(c+c0)
   endif
   ! previous "emax" is replaced by self%e_sc (in namelist: E_storage_capacity)
   ! previous "er" is replaced by sigma_e (see Hense&Beckmann(2006))
   ! previous "fq" is replaced by sigma_q (see Hense&Beckmann(2006))
   ep       = max(0.0_rk,tanh(self%scale-self%scale*e/self%e_sc))
   sigma_e  = tanh(           self%scale*e/self%e_sc)
   fqc      = max(0.0_rk,tanh(self%scale-self%scale*q/self%Qc))
   sigma_q  = tanh(           self%scale*q/self%Qc)

   ! specific rates 
   light_capture = f_pH * self%lightcapture_factor * sigma_l * ep
   uptake = f_pH * self%uptake_factor * omega * sigma_n * sigma_e * fqc
   growth = f_pH * self%growth_factor * omega * sigma_e * sigma_q
   ! fixation_factor is the energy consumption factor for fixation = 3 for heterocysts
   fixation = 0.0_rk
   if (self%n_fixation) then
     fixation = f_pH * omega * sigma_e
     growth = self%fcy2 * growth
   end if

   c_flux=0.0_rk
   g_flux=0.0_rk
   s_flux=0.0_rk
   if (self%lifecycling) call self%calculate_lifecycle_flux(e,q,c_flux,s_flux,g_flux)

   ! Set temporal derivatives into next lifecycle
   _SET_DD_SYM_(self%id_c,self%id_nextc, c*c_flux)
   _SET_DD_SYM_(self%id_s,self%id_nexts, s*s_flux)
   _SET_DD_SYM_(self%id_g,self%id_nextg, g*g_flux)
   ! Set temporal derivatives
   _SET_DD_SYM_(self%id_c,self%id_detritus,c*self%mort)
   _SET_DD_SYM_(self%id_s,self%id_detritus,s*self%mort)
   _SET_DD_(self%id_g,self%id_g,g*self%mort)
   _SET_PP_(self%id_g,self%id_g,c*light_capture)

   _SET_DD_SYM_(self%id_s,self%id_c,c*growth)
   _SET_DD_SYM_(self%id_g,self%id_s,c*uptake)
   _SET_DD_(self%id_g,self%id_g,c*growth)
   
   ni = max(ni, self%minimum_nitrate)

   ! external nutrients
   _SET_DD_(self%id_nitrate,self%id_nitrate,c*uptake*ni/(ni+am))
   if (self%use_ammonium) then
     _SET_DD_(self%id_ammonium,self%id_ammonium,c*uptake*am/(ni+am))
   end if
   ! fixation
   _SET_DD_SYM_(self%id_g,self%id_c,self%growth_factor*c*fixation)
   _SET_DD_(self%id_g,self%id_g,c*self%m*fixation)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_diagcflux,c*c_flux*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_diaggflux,g*g_flux*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_diagsflux,s*s_flux*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_gr,growth*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_up,uptake*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_lc,light_capture*86400.0_rk)
   if (self%n_fixation) then
     _SET_DIAGNOSTIC_(self%id_nf,self%growth_factor*fixation*86400.0_rk)
     _SET_DIAGNOSTIC_(self%id_nf_mean,self%growth_factor*fixation*86400.0_rk)
   end if
   if (self%use_ph) then
     _SET_DIAGNOSTIC_(self%id_f_pH,f_pH)
   end if
   _SET_DIAGNOSTIC_(self%id_primprod,c*(growth + self%growth_factor*fixation)*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_netprimprod,c*(growth + self%growth_factor*fixation - self%mort - c_flux)*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_nfix,c*(self%growth_factor*fixation)*secs_pr_day)  

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do_ppdd


   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_uhh_clcbase), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
   real(rk) :: c,g,s

   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_c,c)                ! biomass
   _GET_(self%id_g,g)                ! energy
   _GET_(self%id_s,s)                ! quota

   ! in bio_clc code, self%depo was given in units 1/s per mmol-N/m3 biomass,
   ! here, self%depo is of units m/s per mmol-N/m3 biomass. The results should
   ! be comparable, since the bottom layer height in the calibrated setup
   ! is appr. 1.03 m.
   _ADD_BOTTOM_FLUX_(self%id_c,-self%depo*c*c)
   _ADD_BOTTOM_FLUX_(self%id_s,-self%depo*c*s)
   _ADD_BOTTOM_FLUX_(self%id_g,-self%depo*c*g)

   ! Export diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sediment_depo,(self%depo*c*c)*86400.0_rk)
   
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

   subroutine do_bottom_ppdd(self,_ARGUMENTS_DO_BOTTOM_PPDD_)
   class (type_uhh_clcbase), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_PPDD_
   real(rk) :: c,g,s

   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_c,c)                ! biomass
   _GET_(self%id_g,g)                ! energy
   _GET_(self%id_s,s)                ! quota

   _ADD_BOTTOM_FLUX_(self%id_c,self%depo*c*c)
   _ADD_BOTTOM_FLUX_(self%id_s,self%depo*c*s)
   _ADD_BOTTOM_FLUX_(self%id_g,self%depo*c*g)

   ! Export diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sediment_depo,(self%depo*c*c)*86400.0_rk)

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom_ppdd


   end module fabm_uhh_clcbase


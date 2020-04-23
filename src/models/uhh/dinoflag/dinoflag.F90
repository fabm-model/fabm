#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_dinoflag - dinoflagellates lifestage model
!
! !INTERFACE:
   module fabm_uhh_dinoflag
!
! !DESCRIPTION:
!
! The dinoflagellates lifestage model
! as published in Warns et al. (2012, JPR)
!
! !USES:
   use fabm_types
   use fabm_expressions

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uhh_dinoflag
!     Variable identifiers
      type (type_state_variable_id)        :: id_veg,id_gam,id_res,id_ger
      type (type_state_variable_id)        :: id_nitrate,id_ammonium
      type (type_state_variable_id)        :: id_phosphate,id_detritus,id_oxygen
      type (type_dependency_id)            :: id_par,id_temp, id_resdep
      type (type_diagnostic_variable_id)   :: id_sl
      type (type_horizontal_diagnostic_variable_id) :: id_tmat
      type (type_horizontal_diagnostic_variable_id) :: id_meanres_diag
      type (type_global_dependency_id)     :: id_doy
      type (type_bottom_state_variable_id) :: bsv_rsum,bsv_rmax
      type (type_surface_state_variable_id):: ssv_meanllim
      type (type_horizontal_dependency_id) :: bsv_meanllim, id_meanres
      type (type_diagnostic_variable_id)   :: id_gr, id_tvg, id_tgr, id_trg, id_tgv, id_llim, id_nlim, id_tlim, id_glim, id_veg_gam_rat 
 
!     Model parameters
      real(rk) :: mumax
      real(rk) :: vmort
      real(rk) :: rkn
      real(rk) :: rdepo
      real(rk) :: alpha_dfl
      real(rk) :: w_veg
      real(rk) :: w_gam
      real(rk) :: w_res
      real(rk) :: w_ger
      real(rk) :: rkc
      real(rk) :: gr_crit
      real(rk) :: l_crit
      real(rk) :: tmat_crit 
      real(rk) :: mdt
      real(rk) :: trate_veg_gam
      real(rk) :: trate_gam_res
      real(rk) :: trate_res_ger
      real(rk) :: trate_ger_veg
      real(rk) :: qg
      real(rk) :: sr
      real(rk) :: s2
      real(rk) :: s3
      real(rk) :: doy_max
      real(rk) :: tmat_initial
      real(rk) :: minimum_nitrate
      logical  :: use_doy_for_maturation
      logical  :: use_phosphate
      logical  :: use_ammonium
      logical  :: use_oxygen

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: check_surface_state
      procedure :: check_bottom_state
      procedure :: get_light_extinction
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------

   contains


   !> Initialise the dinoflag model
   !
   !>  Here, the uhh_dinoflag namelist is read and variables exported
   !!  by the model are registered with FABM.
   subroutine initialize(self,configunit)
   class (type_uhh_dinoflag), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

   real(rk)     :: tmat_initial=170.
   real(rk)     :: mumax= 0.4
   real(rk)     :: vmort=0.03
   real(rk)     :: rkn=0.5
   real(rk)     :: rdepo=5.6e-6
   real(rk)     :: alpha_dfl=0.18
   real(rk)     :: w_veg=0.0
   real(rk)     :: w_gam=0.0
   real(rk)     :: w_res=-8.0
   real(rk)     :: w_ger=8.0
   real(rk)     :: rkc=0.07  
   real(rk)     :: sr=0.0625 
   real(rk)     :: s2=6.625
   real(rk)     :: s3=8.625 
   real(rk)     :: gr_crit=3.75e-7
   real(rk)     :: l_crit = 0.64
   real(rk)     :: tmat_crit = 200.0
   real(rk)     :: mdt=1800.
   real(rk)     :: trate_veg_gam     = 0.015
   real(rk)     :: trate_gam_res     = 0.06
   real(rk)     :: trate_res_ger     = 0.05
   real(rk)     :: trate_ger_veg     = 0.3
   real(rk)     :: qg                = 0.58
   real(rk)     :: doy_max=100.
   real(rk)     :: minimum_nitrate=0.0_rk
   logical      :: use_doy_for_maturation=.false.

   character(len=64)         :: phosphate_variable
   character(len=64)         :: ammonium_variable
   character(len=64)         :: nitrate_variable
   character(len=64)         :: detritus_variable
   character(len=64)         :: oxygen_variable

  namelist /uhh_dinoflag/ mumax,vmort,     &
                      rkn,rdepo, trate_veg_gam,trate_gam_res,  &
                      trate_res_ger,trate_ger_veg,gr_crit,     &
                      w_veg,w_gam,w_res,w_ger,rkc,alpha_dfl,   &
                      l_crit,tmat_crit, mdt,sr,s2,s3,    &
                      doy_max,use_doy_for_maturation,          &
                      ammonium_variable,nitrate_variable,      &
                      phosphate_variable, detritus_variable,   &
                      oxygen_variable,minimum_nitrate,         &
                      tmat_initial

   nitrate_variable = 'uhh_ergom_split_base_nit'
   ammonium_variable = 'uhh_ergom_split_base_amm'
   phosphate_variable = 'uhh_ergom_split_base_pho'
   detritus_variable = 'uhh_ergom_split_base_det'
   oxygen_variable = 'uhh_ergom_split_base_oxy'

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_dinoflag,err=99)

   ! set dependency switches
   self%use_phosphate = phosphate_variable /= ''
   self%use_ammonium  = ammonium_variable /= ''
   self%use_oxygen    = oxygen_variable /= ''

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%mumax, 'mumax', default=mumax,scale_factor=one_pr_day)
   call self%get_parameter(self%vmort, 'vmort', default=vmort,scale_factor=one_pr_day)
   call self%get_parameter(self%rkn, 'rkn', default=rkn)
   call self%get_parameter(self%rdepo, 'rdepo', default=rdepo,scale_factor=one_pr_day)
   call self%get_parameter(self%alpha_dfl, 'alpha_dfl', default=alpha_dfl,scale_factor=one_pr_day)
   call self%get_parameter(self%w_veg, 'w_veg', default=w_veg,scale_factor=one_pr_day)
   call self%get_parameter(self%w_gam, 'w_gam', default=w_gam,scale_factor=one_pr_day)
   call self%get_parameter(self%w_res, 'w_res', default=w_res,scale_factor=one_pr_day)
   call self%get_parameter(self%w_ger, 'w_ger', default=w_ger,scale_factor=one_pr_day)
   call self%get_parameter(self%rkc, 'rkc', default=rkc)
   call self%get_parameter(self%gr_crit, 'gr_crit', default=gr_crit)
   call self%get_parameter(self%l_crit,  'l_crit', default=l_crit)
   call self%get_parameter(self%tmat_crit, 'tmat_crit', default=tmat_crit)
   call self%get_parameter(self%tmat_initial, 'tmat_initial', default=tmat_initial)
   call self%get_parameter(self%mdt, 'mdt', default=mdt)
   call self%get_parameter(self%sr,  'sr',  default=sr)
   call self%get_parameter(self%s2,    's2',    default=s2)
   call self%get_parameter(self%s3,    's3',    default=s3)
   call self%get_parameter(self%minimum_nitrate, 'minimum_nitrate',  default=minimum_nitrate)
   call self%get_parameter(self%trate_veg_gam, 'trate_veg_gam', default=trate_veg_gam,scale_factor=one_pr_day)
   call self%get_parameter(self%trate_gam_res, 'trate_gam_res', default=trate_gam_res,scale_factor=one_pr_day)
   call self%get_parameter(self%trate_res_ger, 'trate_res_ger', default=trate_res_ger,scale_factor=one_pr_day)
   call self%get_parameter(self%trate_ger_veg, 'trate_ger_veg', default=trate_ger_veg,scale_factor=one_pr_day)
   call self%get_parameter(self%qg, 'qg', default=qg)
   call self%get_parameter(self%use_doy_for_maturation, 'use_doy_for_maturation', default=use_doy_for_maturation)
   call self%get_parameter(self%doy_max, 'doy_max', default=doy_max)

   ! Register state variables
   call self%register_state_variable(self%id_veg,'veg', &
         'mmol n/m**3','dinoflagellates vegetative biomass', &
         minimum=0.0e-7_rk,vertical_movement=self%w_veg)

   call self%register_state_variable(self%id_gam,'gam', &
         'mmol n/m**3','dinoflagellates gametes biomass',  &
         minimum=0.0e-7_rk,vertical_movement=self%w_gam)

   call self%register_state_variable(self%id_res,'res', &
         'mmol n/m**3','dinoflagellates resting biomass',  &
         minimum=0.0e-7_rk,vertical_movement=self%w_res)

   call self%register_state_variable(self%id_ger,'ger', &
         'mmol n/m**3','dinoflagellates germinating biomass',  &
         minimum=0.0e-7_rk,vertical_movement=self%w_ger)

   call self%register_bottom_state_variable(self%bsv_rsum,'rsum', &
         'mmol n/m**3','time-integrated maximum bottom biomass', &
         initial_value=3.0_rk*tmat_initial)

   call self%register_bottom_state_variable(self%bsv_rmax,'rmax', &
         'mmol n/m**3','maximum bottom biomass', initial_value=3.0_rk)

   call self%register_surface_state_variable(self%ssv_meanllim,'meanllim', &
         '1/1','mean light limitation', initial_value=0.5_rk)

   ! Register dependencies on external standard variables
   if (self%use_ammonium) &
     call self%register_state_dependency(self%id_ammonium, 'ammonium_target', 'mmol/m**3','ammonium source')
   call self%register_state_dependency(self%id_nitrate, 'nitrate_target', 'mmol/m**3','nitrate source')
   if (self%use_phosphate) &
     call self%register_state_dependency(self%id_phosphate, 'phosphate_target',  'mmol/m**3','phosphate source')
     
   call self%register_horizontal_dependency(self%bsv_meanllim, 'meanllim', '1/1','bottom meanllim')

   ! Register vertical mean of resting cells
   !call self%register_dependency(self%id_resdep, 'resdep', 'mmol/m**3', 'resting biomass')
   !call self%request_coupling(self%id_resdep, 'res')
   call self%register_dependency(self%id_meanres, vertical_mean(self%id_res))
   call self%register_diagnostic_variable(self%id_meanres_diag, 'res_vertmean','mmol n/m**3','vertical mean resting biomass', source=source_do_bottom)
   
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
!   call self%register_diagnostic_variable(self%id_sl,'sigma_l','', &
!      'light limitation factor', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tmat,'tmat','d', &
      'maturation time', output=output_instantaneous, source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_gr,'gr_veg','1/d', &
      'relative growth rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tvg,'tau_veg_gam','1/d', &
      'transition rate from vegetatives to gametes', output=output_instantaneous)     
   call self%register_diagnostic_variable(self%id_tgr,'tau_gam_res','1/d', &
      'transition rate from gametes to cyst', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_trg,'tau_res_ger','1/d', &
      'transition rate from cyst to germinates', output=output_instantaneous)  
   call self%register_diagnostic_variable(self%id_tgv,'tau_ger_veg','1/d', &
      'transition rate from germinates to vegetatives', output=output_instantaneous)        
   call self%register_diagnostic_variable(self%id_llim,'light_limitation','', &
      'light limitation', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_nlim,'nutrient_limitation','', &
      'nutrient limitation', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tlim,'temp_limitation','', &
      'temperature limitation', output=output_instantaneous)      
   call self%register_diagnostic_variable(self%id_glim,'growth_limitation','', &
      'growth limitation', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_veg_gam_rat,'veg_gam_ratio','', &
      'ratio of gam to veg', output=output_instantaneous)       
      
   return

99 call self%fatal_error('fabm_uhh_dinoflag','Error reading namelist uhh_dinoflag')

   end subroutine initialize



   !> Right hand sides of dinoflag model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

   real(rk) :: ni,am,d,po,par,temp,nut
   real(rk) :: veg,gam,res,ger  
   real(rk) :: tlim,llim,nlim,glim,gr_veg
   real(rk) :: veg_gam_rat
   real(rk) :: tau_res_ger,tau_ger_veg,tau_veg_gam,tau_gam_res
   real(rk) :: doy
   real(rk) :: meanres
   real(rk) :: tmat,rsum,rmax
   logical  :: res_ger_maturation

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_veg,veg)       ! vegetatives biomass
   _GET_(self%id_res,res)       ! resting cells biomass
   _GET_(self%id_ger,ger)       ! germinates biomass
   _GET_(self%id_gam,gam)       ! gametes biomass
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

   ! get and calculate maturation time
   _GET_HORIZONTAL_(self%bsv_rmax,rmax)  ! maximum bottom biomass
   _GET_HORIZONTAL_(self%bsv_rsum,rsum)  ! integrated maximum biomass

   if(rmax .eq. 0) then 
     tmat = 0_rk
   else
     tmat = rsum/rmax   
   end if
   
   ! some diagnostics
   
   if(veg .eq. 0) then 
     veg_gam_rat = 0_rk
   else
     veg_gam_rat = gam/veg  
   end if

   ! vegetatives growth
   llim=1.0_rk - exp( -(self%alpha_dfl * par) / self%mumax)
   nlim=nut/(nut+self%rkn)
   tlim=0.5_rk *(tanh(0.8_rk * (temp - 2.3_rk)) - (tanh(3.0_rk *(temp - 8.5_rk)) ))
   glim=0.27_rk *  (2.7_rk -tanh(100.0_rk * (veg_gam_rat - self%qg)))
   gr_veg= self%mumax*llim*nlim*tlim*glim

   ! lifestage fluxes
   tau_res_ger=0.0_rk
   tau_ger_veg=0.0_rk
   tau_veg_gam=0.0_rk
   tau_gam_res=0.0_rk
   !  so far, replaced integrated variables by local conditions
   
   ! ger -> veg
   tau_ger_veg=self%trate_ger_veg * (1.0_rk-exp(-0.5_rk*par))
   ! veg -> gam
   tau_veg_gam=self%trate_veg_gam*(1.0_rk + 4.5_rk*( 0.5_rk*(1.0_rk+tanh((temp-4.5_rk)/0.5_rk)) ) )
   ! gam -> res
   if(veg_gam_rat .gt. self%qg) &
     tau_gam_res=self%trate_gam_res*0.5_rk*(1.0_rk+tanh((temp-4.5_rk)/0.5_rk))
                
   ! Set temporal derivatives
   _ADD_SOURCE_(self%id_veg,veg*(gr_veg - self%vmort) - veg*tau_veg_gam + ger*tau_ger_veg)
   _ADD_SOURCE_(self%id_gam,-gam*self%vmort - gam*tau_gam_res + veg*tau_veg_gam)
   _ADD_SOURCE_(self%id_ger,-ger*self%vmort - ger*tau_ger_veg) 
   _ADD_SOURCE_(self%id_res,gam*tau_gam_res)
   ni = max(ni, self%minimum_nitrate)
   
   ! external nutrients
   _ADD_SOURCE_(self%id_nitrate,-veg*gr_veg * ni/(ni+am))
   if (self%use_ammonium) then
     _ADD_SOURCE_(self%id_ammonium,-veg*gr_veg * am/(ni+am))
   end if
   _ADD_SOURCE_(self%id_detritus,(veg+gam+ger)*self%vmort)
   if (self%use_phosphate) then
     _ADD_SOURCE_(self%id_phosphate, -self%sr *veg*gr_veg)
   end if
   ! add oxygen dynamics
    _ADD_SOURCE_(self%id_oxygen, (self%s2 *am/(am+ni) + self%s3* ni/(am+ni))*veg*gr_veg)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_gr,gr_veg*secs_pr_day)   
   _SET_DIAGNOSTIC_(self%id_tvg,tau_veg_gam*secs_pr_day)   
   _SET_DIAGNOSTIC_(self%id_tgr,tau_gam_res*secs_pr_day)   
   _SET_DIAGNOSTIC_(self%id_trg,tau_res_ger*secs_pr_day)   
   _SET_DIAGNOSTIC_(self%id_tgv,tau_ger_veg*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_llim,llim)
   _SET_DIAGNOSTIC_(self%id_nlim,nlim)
   _SET_DIAGNOSTIC_(self%id_tlim,tlim)
   _SET_DIAGNOSTIC_(self%id_glim,glim)   
   _SET_DIAGNOSTIC_(self%id_veg_gam_rat,veg_gam_rat)
   _GET_HORIZONTAL_(self%id_meanres,meanres)
   !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_meanres_diag,meanres) 
   
   !_SET_DIAGNOSTIC_(self%id_sl,sigma_l)
   

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

   

   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
   class (type_uhh_dinoflag), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
   
   real(rk)                     :: veg,gam,res,ger

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_veg,veg)       ! vegetatives biomass
   _GET_(self%id_res,res)       ! resting cells biomass
   _GET_(self%id_ger,ger)       ! germinates biomass
   _GET_(self%id_gam,gam)       ! gametes biomass

   ! Self-shading
   _SET_EXTINCTION_(self%rkc*(veg+res+ger+gam))

   _LOOP_END_

   end subroutine get_light_extinction


   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
   real(rk),save :: res,rmax,rsum, tmat
   real(rk)      :: meanres, temp, tau_res_ger

   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_res,res)                    ! biomass
   _GET_HORIZONTAL_(self%bsv_rmax,rmax)      ! maximum bottom biomass
   _GET_HORIZONTAL_(self%bsv_rsum,rsum)      ! maximum integrated bottom biomass
!   _ADD_BOTTOM_FLUX_(self%id_res,-self%rdepo*res*res)

   _ADD_BOTTOM_SOURCE_(self%bsv_rsum,rmax*one_pr_day)
   
   if(rmax .eq. 0) then
     tmat = 0_rk
   else
     tmat = rsum/rmax
   end if
   tau_res_ger=0.0_rk

   if ( tmat .gt. self%tmat_crit )&
    tau_res_ger=self%trate_res_ger *((0.5_rk *tanh(8._rk  *(temp+0.1_rk)) +0.5_rk) &
                 -(0.5_rk *tanh(0.8_rk *( temp-6.6_rk)) +0.5_rk))

   _ADD_BOTTOM_FLUX_(self%id_ger, res*tau_res_ger)
   _ADD_BOTTOM_FLUX_(self%id_res,-self%rdepo*res*res - res*tau_res_ger)      
   _GET_HORIZONTAL_(self%id_meanres,meanres)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_meanres_diag,meanres) 
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tmat,tmat)

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom



   subroutine check_bottom_state(self,_ARGUMENTS_CHECK_BOTTOM_STATE_)
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_CHECK_BOTTOM_STATE_
   real(rk) :: rmax,res
   real(rk) :: meanllim
   real(rk) :: doy
   real(rk) :: meanres

   _HORIZONTAL_LOOP_BEGIN_
   ! set rmax
   _GET_(self%id_res,res)                    ! bottom biomass
   _GET_HORIZONTAL_(self%id_meanres,meanres) ! vertically averaged biomass
   _GET_HORIZONTAL_(self%bsv_rmax,rmax)      ! previous maximum biomass
   _SET_HORIZONTAL_(self%bsv_rmax, max(meanres,rmax))

   _GET_HORIZONTAL_(self%bsv_meanllim,meanllim) ! mean light limitation
   _GET_GLOBAL_(self%id_doy,doy)       ! day of year

   ! alternative reset condition:
   if ((meanllim .gt. self%l_crit).and.(doy .lt. 215.)) then
   !if (meanllim .gt. self%l_crit) then
     _SET_HORIZONTAL_(self%bsv_rmax, 0.0_rk)
     _SET_HORIZONTAL_(self%bsv_rsum, 0.0_rk)
   end if
   
   _HORIZONTAL_LOOP_END_
     
   end subroutine check_bottom_state



   subroutine check_surface_state(self,_ARGUMENTS_CHECK_SURFACE_STATE_)
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_CHECK_SURFACE_STATE_
   real(rk)            :: par,llim,meanllim,doy
   ! set averaging window (0.5 * original averaging window, because
   ! of dampening, recursive filter)
   real(rk), parameter :: window = 3.5_rk*86400.0_rk
   real(rk), parameter :: one_pr_window = 1.0_rk/window

   _HORIZONTAL_LOOP_BEGIN_
   _GET_HORIZONTAL_(self%ssv_meanllim, meanllim) ! mean light limitation
   _GET_GLOBAL_(self%id_doy,doy)

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)  ! local photosynthetically active radiation

   llim=1.0_rk - exp(-(self%alpha_dfl * par) / self%mumax)
   meanllim = one_pr_window*((window - self%mdt)*meanllim + self%mdt*llim)
   _SET_HORIZONTAL_(self%ssv_meanllim,meanllim)

   _HORIZONTAL_LOOP_END_

   end subroutine check_surface_state


   end module fabm_uhh_dinoflag


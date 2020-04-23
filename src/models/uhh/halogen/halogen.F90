#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_halogen - generic halogen model
!
! !INTERFACE:
   module fabm_uhh_halogen
!
! !DESCRIPTION:
!
! 
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uhh_halogen
!     Variable identifiers
      type (type_state_variable_id)         :: id_halo
      type (type_state_variable_id)         :: id_doc, id_oxygen, id_det
      type (type_diagnostic_variable_id)    :: id_prod, id_loss
      type (type_horizontal_diagnostic_variable_id) :: id_surface_flux
      type (type_dependency_id)             :: id_phyprod,id_phyloss, id_pH
      type (type_dependency_id)             :: id_par,id_temp,id_uv
      type (type_horizontal_dependency_id)  :: id_pp_air, id_wind

!     Model parameters
      real(rk) :: pp_air
      real(rk) :: beta
      real(rk) :: l_con
      real(rk) :: l_rem
      real(rk) :: l_nit
      real(rk) :: l_uv
      real(rk) :: D_ref !one_pr_Dref
      real(rk) :: I_ref !one_pr_Iref
      real(rk) :: k_photo
      real(rk) :: doc_const
      logical  :: use_ph
      logical  :: use_oxygen
      logical  :: use_phyprod
      logical  :: use_phyloss
      logical  :: read_pp_air
      integer  :: halogen_type
      integer  :: CH3I=1,CHBR3=2

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------

   contains


   !> Initialise the halogen model
   !
   !>  Here, the uhh_halogen namelist is read and variables exported
   !!  by the model are registered with FABM.
   subroutine initialize(self,configunit)
   class (type_uhh_halogen), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

   real(rk)     :: pp_air=-99.
   real(rk)     :: beta=0.0_rk   ! nmol/mmolN
   real(rk)     :: l_con=0.0_rk
   real(rk)     :: l_rem=0.0_rk
   real(rk)     :: l_nit=0.0_rk
   real(rk)     :: l_uv =0.0_rk
   real(rk)     :: D_ref=0.13_rk ! mmolN/m3
   real(rk)     :: I_ref=75.0_rk ! W/m2
   real(rk)     :: doc_const=40.0_rk ! mmolC/m3
   real(rk)     :: k_photo=0.0_rk ! m2/mmolCH3I/mmolC/W/s
   logical      :: use_mean_growth=.false.
   integer      :: halogen_type=1

   character(len=64)         :: ph_variable
   character(len=64)         :: phyprod_variable
   character(len=64)         :: phyloss_variable
   character(len=64)         :: oxygen_variable
   character(len=64)         :: pp_air_variable
   character(len=64)         :: detritus_variable
   character(len=64)         :: nutrient_variable
   character(len=64)         :: uv_variable
   character(len=64)         :: doc_variable

  namelist /uhh_halogen/ pp_air, beta,     &
                      l_con, l_rem, l_nit, l_uv, D_ref, I_ref,&
                      pp_air_variable, &   
                      phyprod_variable,phyloss_variable,      &
                      halogen_type, ph_variable,              &
                      k_photo, doc_const,                     &
                      oxygen_variable, detritus_variable

   ph_variable = ''
   oxygen_variable = ''
   pp_air_variable = ''
   phyloss_variable = ''
   phyprod_variable = ''
   detritus_variable = ''
   doc_variable = ''
   nutrient_variable = ''
   uv_variable = 'uhh_uv_rad'

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_halogen,err=99)

   ! set dependency switches
   self%use_ph = ph_variable /= ''
   self%use_oxygen     = oxygen_variable /= ''
   self%read_pp_air    = pp_air_variable /= ''
   self%use_phyprod    = phyprod_variable /= ''
   self%use_phyloss    = phyloss_variable /= ''

   call self%get_parameter(self%pp_air, 'pp_air', default=pp_air)  

   call self%get_parameter(self%halogen_type, 'halogen_type', default=halogen_type)
   call self%get_parameter(self%beta, 'beta', default=beta)

   call self%get_parameter(self%l_con, 'l_con', default=l_con)

   call self%get_parameter(self%l_rem, 'l_rem', default=l_rem)  

   call self%get_parameter(self%l_nit, 'l_nit', default=l_nit)  
   call self%get_parameter(self%l_uv, 'l_uv', default=l_uv)

   call self%get_parameter(self%D_ref, 'D_ref', default=1.0_rk/D_ref) 

   call self%get_parameter(self%I_ref, 'I_ref', default=1.0_rk/I_ref) 

   call self%get_parameter(self%doc_const, 'doc_const', default=doc_const) 

   call self%get_parameter(self%k_photo, 'k_photo', default=k_photo)

   ! Register state variables
   call self%register_state_variable(self%id_halo,'halogen', &
         'nmol/m**3',' concentration')

   ! Register dependencies on external standard variables
   if (self%use_ph) then
     call self%register_dependency(self%id_ph, 'ph_target','','ph diagnostic value')
     call self%request_coupling(self%id_ph, ph_variable)
   end if
   if (self%use_phyprod) then
     call self%register_dependency(self%id_phyprod, 'phyprod_target','mmol-C/m**3/d','primary production')
     call self%request_coupling(self%id_phyprod, phyprod_variable)
   end if
   if (self%use_phyloss) then
     call self%register_dependency(self%id_phyloss, 'phyloss_target','mmol-C/m**3/d','phytoplankton loss')
     call self%request_coupling(self%id_phyprod, phyprod_variable)
   end if
   if (self%use_oxygen) then
     call self%register_state_dependency(self%id_oxygen, 'oxygen_target','mmol-O2/m**3','dissolved oxygen pool')
     call self%request_coupling(self%id_phyprod, phyprod_variable)
   end if
   if (detritus_variable/='') then
     call self%register_state_dependency(self%id_det, 'detritus_target','mmol-N/m**3','detritus pool')
     call self%request_coupling(self%id_det, detritus_variable)
   end if
   if (doc_variable/='') then
     call self%register_state_dependency(self%id_det, 'doc_target','mmol-C/m**3','doc pool')
     call self%request_coupling(self%id_doc, doc_variable)
   end if

   ! Dependency to UV model
   call self%register_dependency(self%id_uv, 'uv_target','W/m2','UV irradiance')
   call self%request_coupling(self%id_uv, uv_variable)

   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_wind,standard_variables%wind_speed)
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   if (self%read_pp_air) call self%register_dependency(self%id_pp_air,pp_air_variable,'ppm','partial pressure in air')

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_prod,'prod','nmol/m**3/s', &
         ' production rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_loss,'loss','nmol/m**3/s', &
         ' loss rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_surface_flux,'flux','nmol/m**2/s', &
         ' surface flux', output=output_instantaneous)

   return
99 call self%fatal_error('fabm_uhh_halogen','Error reading namelist uhh_halogen')

   end subroutine initialize



   !> Right hand sides of halogen model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uhh_halogen), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

   real(rk) :: halo,oxy,ph,I,temp,D
   real(rk) :: phyprod, phyloss
   real(rk) :: photolysis, const_decay, remin_decay, prod
   real(rk) :: photoprod, rad, doc

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_halo,halo)       ! halogen
   if (_AVAILABLE_(self%id_det)) _GET_(self%id_det,D) ! detritus
   if (_AVAILABLE_(self%id_phyprod)) _GET_(self%id_phyprod,phyprod) ! phytoplankton production
   if (_AVAILABLE_(self%id_phyloss)) _GET_(self%id_phyloss,phyloss) ! phytoplankton loss
   if (_AVAILABLE_(self%id_doc)) then
     _GET_(self%id_doc,doc) ! dissolved organic carbon
   else
     doc = self%doc_const
   end if

   if (self%use_oxygen) then
     _GET_(self%id_oxygen,oxy)
   else
     oxy = 200.0_rk
   end if

   if (self%use_ph) then
     _GET_(self%id_ph,ph)
   else
     ph = 8.0_rk
   end if

   ! Retrieve current environmental conditions.
   _GET_(self%id_uv, I)                ! local UV irradiance
   _GET_(self%id_temp,temp)            ! local temperature

   ! production coupled to primary production
   prod = self%beta * phyprod
   ! photoproduction
   photoprod = self%k_photo * I * doc
   
   ! photolysis
   photolysis = self%l_uv * I*self%I_ref * halo
   ! constant decay
   const_decay = self%l_con * halo
   ! decay coupled to remineralisation
   remin_decay = self%l_rem * self%D_ref * D * 1.066_rk**temp * halo
   
   ! Set temporal derivatives
   _ADD_SOURCE_(self%id_halo,prod + photoprod - photolysis - const_decay - remin_decay)

   ! Set diagnostic variables
   _SET_DIAGNOSTIC_(self%id_prod,prod+photoprod)
   _SET_DIAGNOSTIC_(self%id_loss,photolysis + const_decay + remin_decay)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do



   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_uhh_halogen), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
   real(rk)            :: temp,flux,pp_air,halogen, pp
   real(rk)            :: kw, equilibrium_halogen, wind, schmidt_number

   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)     ! local temperature [degC]
   _GET_(self%id_halo, halogen) ! halogen concentration [nmol/m3]
   _GET_HORIZONTAL_(self%id_wind,wind)     ! wind speed [m/s]
   if (self%read_pp_air) then
     _GET_HORIZONTAL_(self%id_pp_air,pp_air) ! partial pressure 
   else
     pp_air = self%pp_air
   end if

   if (self%halogen_type == self%CH3I) then
     ! Moore & Groszko 1999, Wilke & Chang 1955, Wanninkhof 1992
     schmidt_number = ((62.9_rk/52.9_rk)**0.6_rk) * (2004.0_rk - 93.5_rk*temp + 1.39_rk*temp**2)
     equilibrium_halogen = self%pp_air / exp(13.32_rk - 4338.0_rk/(temp + 273.0_rk))

   elseif (self%halogen_type == self%CHBR3) then
     ! Quack & Wallace 2004
     schmidt_number = 4662.8_rk - 319.45_rk*temp + 9.9012_rk*temp**2 - 0.1159*temp**3
     ! Hense & Quack 2007: 75% transfer velocity compared to CO2
     schmidt_number = 1173.33_rk
     ! Moore et al. 1995:
     equilibrium_halogen = self%pp_air / exp(13.16_rk - 4973.0_rk/(temp + 273.0_rk))

   else
     ! invalid halogen type
   end if
   ! transfer velocity (Nightingale et al. 2000)
   kw = (6.16e-7_rk*wind**2 + 9.25e-7_rk*wind) * sqrt(600.0_rk / schmidt_number)

   flux = kw * (equilibrium_halogen - halogen)

   _ADD_SURFACE_FLUX_(self%id_halo,flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_surface_flux,flux)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface


   end module fabm_uhh_halogen

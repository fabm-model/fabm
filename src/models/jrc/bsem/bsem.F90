#include "fabm_driver.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                        !
!                                                                        !
!                                                                        !
!     BS ecosystem  model developed by Temel Oguz et al                  !
!     Incorpored into fabm by JRC marine modelling team (October 2016)   !
!                     													         !
!                                                                        !
!                                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm\_bsem --- BS biogeochemical model
!
! !INTERFACE:
   MODULE jrc_bsem
!
! !DESCRIPTION:
! The biogeochemical model for the Black Sea pelagic ecosystem proposed  by
! Temel Oguz. It includes two phytoplankton types, four zooplankton types,
! PON, ammonium, nitrate, dissolved oxygen and HS. A detailed description
! could be found in Oguz et al., 2002 and successive papers
! Model currency is  mmol N\,m$^{-3}$,
!
! !USES:
   use fabm_types
!
   implicit none
!
   private
!
! !REVISION HISTORY:
!  Author(s):
! 
!
! !PUBLIC_DERIVED_TYPES:
!
   type,extends(type_base_model), public :: type_jrc_bsem
      ! Variable identifiers
      type (type_state_variable_id)        :: id_pl,id_ps,id_dn,id_zs,id_zl,id_zn,id_zg,id_ni,id_am,id_o2,id_hs
      type (type_bottom_state_variable_id) :: id_fl
      type (type_dependency_id)            :: id_par,id_temp,id_salt,id_rho
      type (type_horizontal_dependency_id) :: id_I_0,id_wind,id_taub
      type (type_diagnostic_variable_id)   :: id_dPAR,id_PPR
      type (type_global_dependency_id)     :: id_yearday
      ! Model parameters
      real(rk) :: sfl_am,sfl_ni,pl0,ps0,alpha_l,alpha_s,sigma_l,sigma_s,beta_l,beta_s,kb,ka_l,ka_s,kn_l,kn_s,mps,mpl,Q10
      real(rk) :: phyZ,g_zs,g_zl,g_zn,g_zg,k_zs,k_zl,k_zn,k_zg,mzs,mzl0,mzn,mnzg,mpzg0,mu_zs,mu_zl,mu_zn,mu_zg
      real(rk) :: epsilon_n0,R0,w_dn,r_n0,r_a0,r_s,r_o,r_u,s1,s2,s4,lds,lsd,tau_crit,lsa,bsa,pvel_c
      real(rk) :: temp_bio_high,temp_bio_low,par_lim,temp_zg_lim,ox_zoop_grazing,ox_min,ox_lim_exc,ox_suboxic,zl_pred_lost

      logical  :: fluff=.false.

      contains

      procedure :: initialize
      procedure :: do
      procedure :: get_light_extinction
      procedure :: do_surface
      procedure :: do_bottom
   end type

   real(rk),parameter           :: secs_per_day=86400._rk
!EOP
!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bsem model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, parameter values are read and variables exported by the model are registered with FABM
!
! !INPUT PARAMETERS:
   class(type_jrc_bsem), intent(inout),target :: self
   integer,              intent(in)           :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%get_parameter(self%fluff,'fluff',                     '-',           'include fluff layer', default=.false.)
   call self%get_parameter(self%sfl_ni,'sfl_ni',                   'mmol N/m2/d', 'constant surface nitrate flux', default=0.0_rk)
   call self%get_parameter(self%sfl_am,'sfl_am',                   'mmol N/m3',   'mmolminimum phyto small concentration??',  default=0.0_rk)
   call self%get_parameter(self%pl0,'pl0',                         'mmol N/m3',   'minimum phyto large concentration', default=0.0225_rk)
   call self%get_parameter(self%ps0,'ps0',                         'mmol N/m3',   'minimum phyto small concentration', default=0.0225_rk)
   call self%get_parameter(self%alpha_l,'alpha_l',                 '-',           'initial slope of P-I curve for large phyto', default=0.8_rk)
   call self%get_parameter(self%alpha_s,'alpha_s',                 '-',           'initial slope of P-I curve for small phyto', default=0.35_rk)
   call self%get_parameter(self%sigma_l,'sigma_l',                 '1/d',         'maximum growth rate large phyto', default=1.2_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%sigma_s,'sigma_s',                 '1/d',         'maximum growth rate small phyto', default=1.0_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%beta_l,'beta_l',                   '-',           'initial slope of P-I curve for large phyto', default=0.0015_rk)
   call self%get_parameter(self%beta_s,'beta_s',                   '-',           'initial slope of P-I curve for small phyto', default=0.35_rk)
   call self%get_parameter(self%kb,'kb',                           '-',           'shelf shading attenuation', default=0.01_rk)
   call self%get_parameter(self%ka_l,'ka_l',                       '-',           'half-saturation for amm uptake large phyto', default=0.3_rk)
   call self%get_parameter(self%ka_s,'ka_s',                       '-',           'half-saturation for amm uptake small phyto', default=0.2_rk)
   call self%get_parameter(self%kn_l,'kn_l',                       '-',           'half-saturation for ni uptake large phyto', default=0.5_rk)
   call self%get_parameter(self%kn_s,'kn_s',                       '-',           'half-saturation for ni uptake small phyto', default=0.3_rk)
   call self%get_parameter(self%mpl,'mpl',                         '1/d',         'mortality rate for large phyto', default=0.005_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%mps,'mps',                         '1/d',         'mortality rate for small phyto', default=0.006_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%Q10,'Q10',                         '-',           'temperature control of large phytoplankton growth', default=2.0_rk)
   call self%get_parameter(self%phyZ,'phyZ',                       '-',           'Phytoplankton assimilation efficiency',  default=0.7_rk)
   call self%get_parameter(self%g_zs,'g_zs',                       '1/d',         'maximum grazing rate of zoo small', default=0.8_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%g_zl,'g_zl',                       '1/d',         'maximum grazing rate of zoo large', default=0.5_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%g_zn,'g_zn',                       '1/d',         'maximum grazing rate of noctiluca', default=0.5_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%g_zg,'g_zg',                       '1/d',         'maximum grazing rate of gelatinous zoo', default=0.15_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%k_zs,'k_zs',                       '-',           'half-sat cte for zoo small grazing', default=0.4_rk)
   call self%get_parameter(self%k_zl,'k_zl',                       '-',           'half-sat cte for zoo large grazing', default=0.5_rk)
   call self%get_parameter(self%k_zn,'k_zn',                       '-',           'half-sat cte for zoo noctiuca grazing', default=0.4_rk)
   call self%get_parameter(self%k_zg,'kn_zg',                      '-',           'half-sat cte for zoo gelatinous grazing', default=0.25_rk)
   call self%get_parameter(self%mzs,'mzs',                         '1/d',         'mortality rate of zoo small', default=0.1_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%mzl0,'mzl0',                       '1/d',         'default mortality rate of zoo large', default=0.25_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%mzn,'mzn',                         '1/d',         'mortality rate of noctiluca', default=0.15_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%mnzg,'mnzg',                       '1/d',         'mortality rate of gelatinous', default=0.02_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%mpzg0,'mpzg0',                     '1/d',         'default predation mortality rate of zoo gelatinous', default=0.1_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%mu_zs,'mu_zs',                     '1/d',         'excretion rate of zoo small', default=0.06_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%mu_zl,'mu_zl',                     '1/d',         'excretion rate of zoo large', default=0.05_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%mu_zn,'mu_zn',                     '1/d',         'excretion rate of zoo noctiluca', default=0.06_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%mu_zg,'mu_zg',                     '1/d',         'excretion rate of zoo gelatinous', default=0.08_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%epsilon_n0,'epsilon_n0',           '1/d',         'default remineralisation rate of dn', default=0.05_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%R0,'R0',                           '-',           'half-sat value for dn remineralisation', default=150._rk)
   call self%get_parameter(self%w_dn,'w_dn',                       'm/d',         'Detritus sedimentation rate', default=5._rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%r_n0,'r_n0',                       '1/d',         'default nitrification rate', default=0.1_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%r_a0,'r_a0',                       '1/d',         'Amm oxid rate by nit', default=0.1_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%r_s,'r_s',                         '1/d',         'HS oxid rate by nit', default=0.1_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%r_o,'r_o',                         '1/d',         'HS oxid rate by oxyg', default=0.1_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%r_u,'r_u',                         '1/d',         'HS oxid rate by other procc', default=0.1_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%s1,'s1',                           '-',           'reduced nitrate/oxidized detritus', default=5.3_rk)
   call self%get_parameter(self%s2,'s2',                           '-',           'oxygen production/recycled nitrogen', default=6.625_rk)
   call self%get_parameter(self%s4,'s4',                           '-',           'nitrification', default=6.625_rk)
   call self%get_parameter(self%lds,'lds',                         'm/d',         'rate of detritus sinking into sediment', default=3.5_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%lsd,'lsd',                         '1/d',         'rate of sediment resuspension', default=25._rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%tau_crit,'tau_crit',               'N /m2',       'critical bottom stress',  default=0.07_rk)
   call self%get_parameter(self%lsa,'lsa',                         '1/d',         'rate of sediment mineralisation', default=0.001_rk, scale_factor=1._rk/secs_per_day)
   call self%get_parameter(self%bsa,'bsa',                         '1/deg C',     'temperature control of sediment mineralisation',  default=0.15_rk)
   call self%get_parameter(self%pvel_c,'pvel_c',                   '-',           'A constant for adjusting oxygen flux', default=3.5_rk)
   call self%get_parameter(self%temp_bio_high,'temp_bio_high',     'deg C',       'Upper temperature control of biological reactions',  default=20.0_rk)
   call self%get_parameter(self%temp_bio_low,'temp_bio_low',       'deg C',       'Lower temperature control of biological reactions',  default=5.0_rk)
   call self%get_parameter(self%par_lim,'par_lim',                 'W/m**2',      'PAR limitation for photo inhibition',  default=75.0_rk)
   call self%get_parameter(self%temp_zg_lim,'temp_zg_lim',         'deg C',       'Gelationous zooplankton critical temperature' ,  default=16.0_rk)
   call self%get_parameter(self%ox_zoop_grazing,'ox_zoop_grazing', 'mmol O2/m3',  'Oxygen control for zooplankton grazing',  default=200.0_rk)
   call self%get_parameter(self%ox_min,'ox_min',                   'mmol O2/m3',  'Oxygen minimum allowed value',  default=0.05_rk)
   call self%get_parameter(self%ox_lim_exc,'ox_lim_exc',           'mmol O2/m3',  'Oxygen control for zooplankton excretion',  default=300.0_rk)
   call self%get_parameter(self%ox_suboxic,'ox_suboxic',           'mmol O2/m3',  'Oxygen control for suboxic conditions',  default=10.0_rk)
   call self%get_parameter(self%zl_pred_lost,'zl_pred_lost',       '-',           'Fraction of Zl lost by predation',  default=0.5_rk)
 
   ! Register state variables
   call self%register_state_variable(self%id_pl,'pl','mmol N/m**3','PhytoLarge',minimum=0.0_rk)
   call self%register_state_variable(self%id_ps,'ps','mmol N/m**3','PhytoSmall',minimum=0.0_rk)
   call self%register_state_variable(self%id_zs,'zs','mmol N/m**3','MicroZoo',minimum=0.0_rk)
   call self%register_state_variable(self%id_zl,'zl','mmol N/m**3','MesoZoo',minimum=0.0_rk)
   call self%register_state_variable(self%id_zn,'zn','mmol N/m**3','Noctiluca',minimum=0.0_rk)
   call self%register_state_variable(self%id_zg,'zg','mmol N/m**3','Mnemiopsis',minimum=0.0_rk)
   call self%register_state_variable(self%id_dn,'dn','mmol N/m**3','detritus',minimum=0.0_rk,vertical_movement=self%w_dn)
   call self%register_state_variable(self%id_am,'am','mmol N/m**3','ammonium',minimum=0.0_rk)
   call self%register_state_variable(self%id_ni,'ni','mmol N/m**3','nitrate',minimum=0.0_rk)
   call self%register_state_variable(self%id_o2,'o2','mmol O2/m**3','oxygen')
   call self%register_state_variable(self%id_hs,'hs','mmol HS/m**3','hydrogen sulfide')
   if (self%fluff) call self%register_state_variable(self%id_fl,'fl','mmol N/m**2','flf',minimum=0.0_rk)

   ! Let all nitrogen-based state variables contribute to total nitrogen to enable mass conservation checking.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_pl)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_ps)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_zs)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_zl)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_zn)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_zg)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_dn)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_am)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_ni)
   if (self%fluff) call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_fl)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2',       'photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_PPR, 'PPR','mmol N/m**3/d','gross primary production')

   ! Register environmental dependencies
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_rho,standard_variables%density)
   call self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_wind,standard_variables%wind_speed)
   call self%register_dependency(self%id_yearday,standard_variables%number_of_days_since_start_of_the_year)

   if (self%fluff) call self%register_dependency(self%id_taub,standard_variables%bottom_stress)

#if 0
write(*,*) self%sfl_am,self%sfl_ni,self%fluff,self%pl0,self%ps0
write(*,*) self%alpha_l,self%alpha_s,self%sigma_l,self%sigma_s
write(*,*) self%beta_l,self%beta_s,self%kb,self%ka_l,self%ka_s,self%kn_l,self%kn_s
write(*,*) self%mpl,self%mps,self%Q10,self%phyZ,self%g_zs,self%g_zl
write(*,*) self%g_zn,self%g_zg,self%k_zs,self%k_zl,self%k_zn,self%k_zg,self%mzs
write(*,*) self%mzl0,self%mzn,self%mnzg,self%mpzg0,self%mu_zs,self%mu_zl
write(*,*) self%mu_zn,self%mu_zg,self%epsilon_n0,self%R0
write(*,*) self%w_dn,self%r_n0,self%r_a0,self%r_s,self%r_o,self%r_u,self%pvel_c
write(*,*) self%s1,self%s2,self%s4,self%lds
write(*,*) self%lsd,self%tau_crit,self%lsa,self%bsa
write(*,*) self%zl_pred_lost,self%ox_suboxic,self%ox_lim_exc,self%temp_bio_low,self%temp_bio_high
write(*,*) self%par_lim,self%temp_zg_lim,self%ox_min,self%ox_zoop_grazing
#endif

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of bsem model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!  The right hand sides of the Oguz et al., biogeochemical model are
!  coded in this soubroutine.
!
! !INPUT PARAMETERS:
   class(type_jrc_bsem), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)        :: pl,ps,zs,zl,zn,zg,am,ni,dn,ox,hs,par,I_0,temp,rho,dens,yearday
   real(rk)        :: beta_l0,beta_s0,ppi_l,ppi_s,temp_lim_l,temp_lim_s,am_lim_s,am_lim_l
   real(rk)        :: nit_lim_l,nit_lim_s,growth_lim_pl,growth_lim_ps,zoo_ox_cont,mort_temp
   real(rk),dimension(6):: prey_vect,b_1,b_2,b_3,b_4,g_1,g_2,g_3,g_4
   real(rk),dimension(4):: pred_vect
   real(rk),dimension(6,4):: a,b,G
   real(rk)        :: temp_cont_zg,temp_cont_zn,mpzg0_loc
   real(rk)        :: epsilon_n,r_n,nitrif,amm_nit_ox,r_a,denitri,oxilimdecom,oxilimexc
   real(rk)        :: HS_ox_oxidation,HS_nit_oxidation,HS_oth_oxidation,mzl,mpzg
  
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial_loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values
   _GET_(self%id_pl,pl) !large phyto
   _GET_(self%id_ps,ps) !small phyto
   _GET_(self%id_zs,zs) !microzoo
   _GET_(self%id_zl,zl) !mesozoo
   _GET_(self%id_zn,zn) !nocticula
   _GET_(self%id_zg,zg) !Mniopsis
   _GET_(self%id_dn,dn) !detritus
   _GET_(self%id_am,am) !ammonium
   _GET_(self%id_ni,ni) !nitrate
   _GET_(self%id_o2,ox) !oxygen
   _GET_(self%id_hs,hs) !hidrogen sulfate

   ! Retrieve current environmental conditions
   _GET_   (self%id_par,par) ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,I_0) ! surface short wave radiation
   _GET_   (self%id_temp,temp) !local water temperature
   _GET_   (self%id_rho,rho) !local water density
   _GET_GLOBAL_ (self%id_yearday,yearday) !decimal day of the year

   ! Convert density into sigma
   dens=rho-1000.0_rk

   ! Temperature control for biological reactions
   if (temp .gt. self%temp_bio_high) then
      temp=self%temp_bio_high
   end if
   if (temp .lt. self%temp_bio_low) then
     temp=self%temp_bio_low
   end if

   ! Light limitation of PP based on light intensity
   if (par .lt. self%par_lim) then
      beta_l0=0.0_rk
      beta_s0=0.0_rk
else
      beta_l0 = self%beta_l
      beta_s0 = self%beta_s
   end if
   ppi_l = tanh(self%alpha_l*par)*exp(-beta_l0*par)
   ppi_s = tanh(self%alpha_s*par)*exp(-beta_s0*par)

   ! Temperature limitation of PP
   
      temp_lim_l=self%Q10**((-temp+12.0_rk)/12.0_rk)
      temp_lim_s=1._rk

   ! Nutrients limitation of PP

   am_lim_l=am/(self%ka_l+am)
   am_lim_s=am/(self%ka_s+am)

   nit_lim_l=(ni/(self%kn_l+ni))*(self%ka_l/(self%ka_l+am))
   nit_lim_s=(ni/(self%kn_s+ni))*(self%ka_s/(self%ka_s+am))

   ! Combined growth limitation
  
   growth_lim_pl=self%sigma_l*ppi_l*temp_lim_l*(am_lim_l+nit_lim_l)
   growth_lim_ps=self%sigma_s*ppi_s*temp_lim_s*(am_lim_s+nit_lim_s)

   ! Food preference for each zooplankton group
   prey_vect=reshape([pl,ps,dn,zs,zl,zn],shape(prey_vect))
   pred_vect=reshape([zs,zl,zn,zg],shape(pred_vect))
   a=reshape([0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.5_rk,0.5_rk,0.0_rk,0.25_rk, &
       0.5_rk,0.5_rk,1.0_rk,0.25_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,0.5_rk,0.3_rk],shape(a))
   b_1=(a(:,1)*prey_vect)/sum(a(:,1)*prey_vect)
   b_2=(a(:,2)*prey_vect)/sum(a(:,2)*prey_vect)
   b_3=(a(:,3)*prey_vect)/sum(a(:,3)*prey_vect)
   b_4=(a(:,4)*prey_vect)/sum(a(:,4)*prey_vect)
   b=reshape([b_1,b_2,b_3,b_4],shape(b))

   ! Temperature control for zooplankton grazing
   ! For noctiluca
   temp_cont_zn=self%Q10**((temp-12.0_rk)/8.0_rk)
  
   ! For gelatinous
   if(temp .gt. self%temp_zg_lim) then
     temp_cont_zg=1.0_rk+(3.5_rk*((temp-16.0_rk)/(1.0_rk+(temp-16.0_rk))))
   else
     temp_cont_zg=1.0_rk
   end if

   ! Oxygen control of zooplankton grazing
   if (ox .gt. self%ox_zoop_grazing) then
      zoo_ox_cont=1.0_rk
   end if
   if (ox .lt. self%ox_zoop_grazing  .and. ox .gt. self%ox_min) then
      zoo_ox_cont=ox/(self%R0+ox)
   end if
   if (ox .lt. self%ox_min) then
!AS      zoo_ox_cont=0.05_rk
      zoo_ox_cont=self%ox_min
   end if

   ! Grazing terms
   g_1=self%g_zs*((b(:,1)*prey_vect)/(self%k_zs+sum(b(:,1)*prey_vect)))
   g_2=self%g_zl*zoo_ox_cont*((b(:,2)*prey_vect)/(self%k_zl+sum(b(:,2)*prey_vect)))
   g_3=self%g_zn*zoo_ox_cont*temp_cont_zn*((b(:,3)*prey_vect)/(self%k_zn+sum(b(:,3)*prey_vect)))
   g_4=self%g_zg*temp_cont_zg*((b(:,4)*prey_vect)/(self%k_zg+sum(b(:,4)*prey_vect)))
   G=reshape([g_1,g_2,g_3,g_4],shape(G))

   ! Zoo Mortality terms
   ! Mesozoo mortality dependent on temp
   mort_temp=1.0_rk-((temp-5.0_rk)/(3.0_rk+(temp-5.0_rk)))
   mzl=self%mzl0*mort_temp

   ! Gelatinous zoop mortality dependent on density
   ! and it's time regulation
 
   ! This is to parametrize Beroa predation on Mnemiopsis

   ! Maximum predation mortality is four times larger towards the end of the year
   if (yearday .lt. 300.0_rk) then
      mpzg0_loc=self%mpzg0
   else
      mpzg0_loc=4.0_rk*self%mpzg0
   end if
   
   ! And then, Beroa predation is only happening on the surface layers
   if (dens .lt. 14.0_rk) then
      mpzg=mpzg0_loc/(self%k_zg+zg)
   else
      mpzg=0.0_rk
   end if

   ! Zooplankton excretion control
   if (ox .lt. self%ox_lim_exc) then
      oxilimexc=0.0_rk
   else
      oxilimexc=1.0_rk
   end if

   ! Decomposition of PON

   if(ox .ge. 250.0_rk) then
     oxilimdecom=1.0_rk
   end if 
   if(ox .lt. 250.0_rk .and. ox .gt. self%ox_min ) then
      oxilimdecom=1.0_rk-(self%R0/(self%R0+ox))
   else
      oxilimdecom=0.25_rk
    end if
   
   if (oxilimdecom .lt. 0.25_rk) then
      oxilimdecom=0.25_rk
   end if
   
   epsilon_n=self%epsilon_n0*oxilimdecom

   ! Ammonium nitrification
   if (ox .gt. self%ox_min ) then
      r_n=self%r_n0*(10.0_rk/(10.0_rk+par))
   else
      r_n=0.0_rk
   end if
   nitrif=r_n*am

   ! Ammonium oxidation by nitrate
   if (ox .lt. self%ox_suboxic .and. ox .gt. self%ox_min) then
      r_a=self%r_a0
   else
      r_a=0.0_rk
   end if
   amm_nit_ox=r_a*am*ni

   ! Decomposition of nitrate  ! this is a lost from the system as nitrate is converted into molecular nitrogen
   if (ox .lt. 20.0_rk ) then
      denitri=0.8_rk*epsilon_n*dn
   else
      denitri=0.0_rk
   end if


   ! Dissolved oxigen and HS
   if (ox .gt. self%ox_suboxic) then
      HS_ox_oxidation=0.5_rk*self%r_o*hs*ox
   else
      HS_ox_oxidation=0.0_rk
   end if
   HS_nit_oxidation=self%r_s*hs*ni
   HS_oth_oxidation=self%r_u*hs

 
!
   ! Concentration control of mniopsis
   if (zg .lt. .00000001_rk) then
      mpzg=0.0_rk
   end if

   ! We start with the differential equations
   _ADD_SOURCE_(self%id_pl,growth_lim_pl*pl-sum(G(1,:)*pred_vect)-self%mpl*pl)
   _ADD_SOURCE_(self%id_ps,growth_lim_ps*ps-sum(G(2,:)*pred_vect)-self%mps*ps)
   _ADD_SOURCE_(self%id_zs,(self%phyZ*sum(G(:,1)*zs))-sum(G(4,:)*pred_vect)-self%mu_zs*zs-self%mzs*zs**2)
   _ADD_SOURCE_(self%id_zl,(self%phyZ*sum(G(:,2)*zl))-sum(G(5,:)*pred_vect)-self%mu_zl*zl-     mzl*zl**2)
   _ADD_SOURCE_(self%id_zn,(self%phyZ*sum(G(:,3)*zn))-sum(G(6,:)*pred_vect)-self%mu_zn*zn-self%mzn*zn**2)
   _ADD_SOURCE_(self%id_zg,(self%phyZ*sum(G(:,4)*zg))-self%mu_zg*zg-self%mnzg*zg-mpzg*zg**2)
   _ADD_SOURCE_(self%id_dn,(1.0_rk-self%phyZ)*(sum(G(:,1)*zs)+sum(G(:,2)*zl)+sum(G(:,3)*zn)+sum(G(:,4)*zg))+self%mpl*pl+self%mps*ps+self%mzs*zs**2+self%zl_pred_lost*(mzl*zl**2)+self%mzn*zn**2+self%mnzg*zg+mpzg*zg**2+(1.0_rk-oxilimexc)*(self%mu_zs*zs+self%mu_zl*zl+self%mu_zn*zn+self%mu_zg*zg)-sum(G(3,:)*pred_vect)-epsilon_n*dn)
   _ADD_SOURCE_(self%id_am,epsilon_n*dn+oxilimexc*(self%mu_zs*zs+self%mu_zl*zl+self%mu_zn*zn+self%mu_zg*zg)-(self%sigma_s*ppi_s*temp_lim_s*am_lim_s*ps)-(self%sigma_l*ppi_l*temp_lim_l*am_lim_l*pl)-nitrif-amm_nit_ox)
   _ADD_SOURCE_(self%id_ni,nitrif-(self%sigma_s*ppi_s*temp_lim_s*nit_lim_s*ps)-(self%sigma_l*ppi_l*temp_lim_l*nit_lim_l*pl)-denitri-0.6_rk*amm_nit_ox-1.34_rk*HS_nit_oxidation)
   _ADD_SOURCE_(self%id_o2,8.1258_rk*(growth_lim_pl*pl+growth_lim_ps*ps)-6.625_rk*(epsilon_n*dn+self%mu_zs*zs+self%mu_zl*zl+self%mu_zn*zn+self%mu_zg*zg)-2.0_rk*nitrif-HS_ox_oxidation)
   _ADD_SOURCE_(self%id_hs,(0.5_rk*epsilon_n*dn)-2.0_rk*HS_ox_oxidation-HS_nit_oxidation-HS_oth_oxidation)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_PPR ,(growth_lim_pl*(pl+self%pl0)+growth_lim_ps*(ps+self%ps0))*secs_per_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

   END subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the bsem model
!
! !INTERFACE:
  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !DESCRIPTION:
!  Here, surface fluxes are computed. For O2 the Wanninkhof (1992) formulation is
!  used. The p-vel is computed as (0.31*wnd**2)*sqrt(660/sc). Being the schmidt
!  number computed as sc=2073.1-125.62*T+3.6276*T**2-0.043219*T**3
!
! !INPUT PARAMETERS:
   class(type_jrc_bsem),intent(in)       ::self
  _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
  real(rk)                 :: temp,wnd,salt,o2,ni,am
  real(rk)                 :: p_vel,sc,flo2
  integer,parameter        :: newflux=1
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_HORIZONTAL_(self%id_wind,wnd)

   _GET_(self%id_o2,o2)
   _GET_(self%id_ni,ni)
   _GET_(self%id_am,am)

   ! Calculation of the surface oxygen flux
!   if (newflux .eq. 1) then
   sc=1450._rk+(1.1_rk*temp-71._rk)*temp
   if (wnd .gt. 13._rk) then
      p_vel = 5.9_rk*(5.9_rk*wnd-49.3_rk)/sqrt(sc)
   else
         if (wnd .lt. 3.6_rk) then
            p_vel = 1.003_rk*wnd/(sc)**(0.66_rk)
         else
            p_vel = 5.9_rk*(2.85_rk*wnd-9.65_rk)/sqrt(sc)
         end if
      end if
      if (p_vel .lt. 0.05_rk) then
         p_vel = 0.05_rk
      end if
      p_vel = p_vel/secs_per_day
      flo2 =self%pvel_c*p_vel*(osat_weiss(temp,salt)-o2)
      _ADD_SURFACE_FLUX_(self%id_o2,flo2)
!   end if

  _ADD_SURFACE_FLUX_(self%id_ni,self%sfl_ni/secs_per_day)
  _ADD_SURFACE_FLUX_(self%id_am,self%sfl_am/secs_per_day)

   _HORIZONTAL_LOOP_END_
   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of benthic_predator model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
!  This routine calculates the benthic sink and source terms, as well as
!  (matching) bottom fluxes for pelagic variables. Both have units mmol/m**2/s.
!  Benthic processes are explained in the description of iow\_ergom\_do
!  subroutine.
!  AT THE MOMENT THIS IS NOT IMPLEMENTED, IT IS EQUAL TO ERGOM MODEL
!
! !INPUT PARAMETERS:
   class (type_jrc_bsem),       intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: fl,amb,nib,pob,deb,oxb,taub,temp
   !logical                   :: fluff,ltaub
   real(rk)                   :: llds,llsd,llsa
   real(rk)                   :: thopnp,thomnp,thomnm,thsum
   real(rk), parameter        :: wo=30.0_rk,wn=0.1_rk,dot2=0.2_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).

   if (self%fluff) then

   _HORIZONTAL_LOOP_BEGIN_
   ! Retrieve current (local) state variable values.
      _GET_(self%id_am,amb)
      _GET_(self%id_dn,deb)
      _GET_(self%id_ni,nib)
      _GET_(self%id_o2,oxb)
      _GET_HORIZONTAL_(self%id_fl,fl)

      _GET_HORIZONTAL_(self%id_taub,taub)
      _GET_(self%id_temp,temp)

      thopnp=th( oxb,wo,0.0_rk,1.0_rk)*yy(wn,nib)
      thomnp=th(-oxb,wo,0.0_rk,1.0_rk)*yy(wn,nib)
      thomnm=th(-oxb,wo,0.0_rk,1.0_rk)*(1.0_rk-yy(wn,nib))
      thsum=thopnp+thomnp+thomnm
      thopnp=thopnp/thsum
      thomnp=thomnp/thsum
      thomnm=thomnm/thsum
      llsa=self%lsa*exp(self%bsa*temp)*(th(oxb,wo,dot2,1.0_rk))

!      ltaub=taub**2*1000.

      if (self%tau_crit .gt. taub) then
         llds=self%lds*(self%tau_crit-taub)/self%tau_crit
      else
         llds=0.0_rk
      end if
      if (self%tau_crit .lt. taub) then
         llsd=self%lsd*(taub-self%tau_crit)/self%tau_crit
      else
         llsd=0.0_rk
      end if

      _ADD_BOTTOM_SOURCE_(self%id_fl,llds*deb-llsd*fl-llsa*fl-th(oxb,wo,0.0_rk,1.0_rk)*llsa*fl)
      _ADD_BOTTOM_FLUX_(self%id_dn,-llds*deb+llsd*fl)
      _ADD_BOTTOM_FLUX_(self%id_am,llsa*fl)
      _ADD_BOTTOM_FLUX_(self%id_ni,self%s1*thomnp*llsa*fl)
      _ADD_BOTTOM_FLUX_(self%id_o2,-(self%s4+self%s2*(thopnp+thomnm))*llsa*fl)

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end if

   end subroutine do_bottom
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
   class (type_jrc_bsem), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                     :: pl,ps,dn
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_pl,pl) ! large phyto
   _GET_(self%id_ps,ps) ! small phyto
   _GET_(self%id_dn,dn) ! detritus

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kb*(self%pl0+self%ps0+pl+ps+dn))

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Weiss formula for the saturation oxygen (osat)
!
! !INTERFACE:
   real(rk) function osat_weiss(t,s)
!
! !DESCRIPTION:
!  Weiss formula for the saturation oxygen (osat) \cite{Weiss1970}
!
! !INPUT PARAMETERS:
  real(rk), intent(in)                 :: t,s
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                 :: tk
   real(rk),parameter       :: aa1=-173.4292_rk
   real(rk),parameter       :: aa2=249.6339_rk
   real(rk),parameter       :: a3=143.3483_rk
   real(rk),parameter       :: a4=-21.8492_rk
   real(rk),parameter       :: b1=-0.033096_rk
   real(rk),parameter       :: b2=0.014259_rk
   real(rk),parameter       :: b3=-0.001700_rk
   real(rk),parameter       :: kelvin=273.16_rk
   real(rk),parameter       :: mol_per_liter=44.661_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   tk=(t+kelvin)*0.01_rk
   osat_weiss=exp(aa1+aa2/tk+a3*log(tk)+a4*tk    &
              +s*(b1+(b2+b3*tk)*tk))*mol_per_liter
   end function osat_weiss
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Step function
!
! !INTERFACE:
   real(rk) function th(x,w,min,max)
!
! !DESCRIPTION:
!  Instead of the
!  heavyside switches used by \cite{Neumannetal2002}, we apply here a smoothed
!  {\it tangens hyperbolicus} transition with prescribed width $x_w$:
!  \begin{equation}\label{theta}
!  \theta (x,x_w,y_{\min},y_{\max})= y_{\min}+(y_{\max}-y_{\min})
!  \frac12\left(1-\tanh \left(\frac{x}{x_w}   \right)      \right).
!  \end{equation}
!
! !INPUT PARAMETERS:
!   type(type_jrc_bsem), INTENT(IN) :: self
   real(rk), intent(in)            :: x,w,min,max
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (w .gt. 1.e-10) then
      th=min+(max-min)*0.5_rk*(1._rk+tanh(x/w))
   else
      if(x .gt. 0.0_rk) then
         th=1.0_rk
      else
         th=0.0_rk
      end if
   end if
   end function th
!EOC

!-----------------------------------------------------------------------
! BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:
   real(rk) function yy(a,x)
!
! !DESCRIPTION:
!  This is a squared Michaelis-Menten type of limiter:
!  \begin{equation}\label{Y}
!  Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
!  \end{equation}
!
! !INPUT PARAMETERS:
  real(rk),intent(in)        :: a,x
!
!EOP
!-----------------------------------------------------------------------
!BOC
   yy=x**2/(a**2+x**2)
   end function yy
!EOC

!-----------------------------------------------------------------------

  end module jrc_bsem

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

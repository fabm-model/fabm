#include "fabm_driver.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                        !
!                                                                        !
!     ERGOM model GOTM version ported to FABM by Gennadi Lessin          !
!     (Marine Systems Institute at Tallinn University of Technology).    !
!     The work performed at EC Joint Research Centre (Ispra).            !
!                                                                        !
!                                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! GOTM implementation of the biogeochemical model ERGOM
!
! The biogeochemical model by
! \cite{Neumannetal2002} consists of $I=10$
! state variables. The nutrient state variables are dissolved
! ammonium, nitrate, and phosphate. Primary production is provided
! by three functional phytoplankton groups: diatoms, flagellates,
! and blue-green algae (cyanobacteria). Diatoms represent larger
! cells which grow fast in nutrient-rich conditions. Flagellates
! represent smaller cells with an advantage at lower nutrients
! concentrations especially during summer conditions. The
! cyanobacteria group is able to fix and utilise atmospheric
! nitrogen and therefore, the model assumes phosphate to be the only
! limiting nutrient for cyanobacteria. Due to the ability of
! nitrogen fixation, cyanobacteria are a nitrogen source for the
! system. A dynamically developing bulk zooplankton variable
! provides grazing pressure on phytoplankton. Dead particles are
! accumulated in a detritus state variable. The detritus is
! mineralised into dissolved ammonium and phosphate during the
! sedimentation process. A certain amount of the detritus reaches
! the bottom, where it is accumulated in the sedimentary detritus.
! Detritus in the sediment is either buried in the sediment,
! mineralised or resuspended into the water column, depending on the
! velocity of near-bottom currents. The development of oxygen in the
! model is coupled to the biogeochemical processes via
! stoichiometric ratios. Oxygen concentration controls processes as
! denitrification and nitrification.
! The basic structure of the model is explained in figure \ref{fig_neumann},
! and a detailed description of the
! model is given in section \ref{sec:bio-gotm-details}.
! \begin{figure}
! \begin{center}
! \scalebox{0.5}{\includegraphics{figures/iow_structure.eps}}
! \caption{Structure of the \cite{Neumannetal2002} model
! with cyanobacteria (cya),
! diatoms (dia), dinoflagellates (fla), detritus (det), zooplankton (zoo),
! ammonium (amm), nitrate (nit) detritus sediment (sed), oxygen (oxy)
! and phosphorus (pho) as the ten
! state variables.
! The concentrations are in mmol N\,m$^{-3}$,
!  mmol N\,m$^{-2}$,  mmol P\,m$^{-3}$ and l O$_2$m$^{-3}$.
! Conservative fluxes are denoted by thin green arrows, non-conservative fluxes
! by bold arrows.
! }\label{fig_neumann}
! \end{center}
! \end{figure}

module gotm_ergom

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_gotm_ergom
      ! Variable identifiers
      type (type_state_variable_id)        :: id_p1,id_p2,id_p3,id_zo,id_de,id_am,id_ni,id_po,id_o2
      type (type_bottom_state_variable_id) :: id_fl
      type (type_dependency_id)            :: id_par,id_temp,id_salt
      type (type_surface_dependency_id)    :: id_I_0,id_wind
      type (type_bottom_dependency_id)     :: id_taub
      type (type_diagnostic_variable_id)   :: id_dPAR,id_GPP,id_NCP,id_PPR,id_NPR

      ! Model parameters
      real(rk) :: sfl_po,sfl_am,sfl_ni,p10,p20,p30,zo0,kc,i_min,r1max,r2max,r3max,alpha1,alpha2,alpha3,lpa,lpd
      real(rk) :: tf,tbg,beta_bg,g1max,g2max,g3max,lza,lzd,iv,topt,lan,oan,beta_an,lda,tda,beta_da
      real(rk) :: pvel,sr,s1,s2,s3,s4,a0,a1,a2,lds,lsd,tau_crit,lsa,bsa,ph1,ph2
      logical  :: fluff
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_gotm_ergom), intent(inout), target :: self
      integer,                 intent(in)            :: configunit

      real(rk) :: w_p1
      real(rk) :: w_p2
      real(rk) :: w_p3
      real(rk) :: w_de
      real(rk), parameter :: secs_pr_day=86400.0_rk

      ! Store parameter values in our own derived type
      ! NB! All rates must be provided in values per day in the configuration file,
      ! and are converted here to values per second by specifying scale_factor=1.0_rk/secs_pr_day
      call self%get_parameter(self%sfl_po, 'sfl_po', 'mmol P/m2/d', 'surface phosphate flux', default=0.0015_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%sfl_ni, 'sfl_ni', 'mmol N/m2/d', 'surface nitrate flux', default=0.09_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%sfl_am, 'sfl_am', 'mmol N/m2/d', 'surface ammonium flux', default=0.07_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%fluff, 'fluff', '', 'use fluff', default=.false.)
      call self%get_parameter(self%p10, 'p10', 'mmol N/m3', 'background concentration of diatoms', default=0.0225_rk)
      call self%get_parameter(self%p20, 'p20', 'mmol N/m3', 'background concentration of flagellates', default=0.0225_rk)
      call self%get_parameter(self%p30, 'p30', 'mmol N/m3', 'background concentration of cyanobacteria', default=0.0225_rk)
      call self%get_parameter(self%zo0, 'zo0', 'mmol N/m3', 'background concentration of zooplankton', default=0.0225_rk)
      call self%get_parameter(self%kc, 'kc', 'm2/mmol N', 'specific light attenuation of phytoplankton and detritus', default=0.03_rk)
      call self%get_parameter(self%i_min, 'i_min', 'W/m2', 'minimum optimal photosynthetically active radiation', default=25.0_rk)
      call self%get_parameter(self%r1max, 'r1max', '1/d', 'maximum growth rate of diatoms', default=2.0_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%r2max, 'r2max', '1/d', 'maximum growth rate of flagellates', default=0.7_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%r3max, 'r3max', '1/d', 'maximum growth rate of cyanobacteria', default=0.5_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%alpha1, 'alpha1', 'mmol N/m3', 'half-saturation nutrient concentration for diatoms', default=1.35_rk)
      call self%get_parameter(self%alpha2, 'alpha2', 'mmol N/m3', 'half-saturation nutrient concentration for flagellates', default=0.675_rk)
      call self%get_parameter(self%alpha3, 'alpha3', 'mmol N/m3', 'half-saturation nutrient concentration for cyanobacteria', default=0.5_rk)
      call self%get_parameter(self%lpa, 'lpa', '1/d', 'phytoplankton exudation', default=0.01_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%lpd, 'lpd', '1/d', 'phytoplankton mortality', default=0.02_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%tf, 'tf', 'deg C', 'half-saturation temperature for flagellates', default=10.0_rk)
      call self%get_parameter(self%tbg, 'tbg', 'deg C', 'reference temperature for cyanobacteria', default=14.0_rk)
      call self%get_parameter(self%beta_bg, 'beta_bg', '1/deg C', 'temperature sensitivity for cyanobacteria', default=1.0_rk)
      call self%get_parameter(self%g1max, 'g1max', '1/d', 'maximum grazing rate on diatoms', default=0.5_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%g2max, 'g2max', '1/d', 'maximum grazing rate on flagellates', default=0.5_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%g3max, 'g3max', '1/d', 'maximum grazing rate on cyanobacteria', default=0.25_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%lza, 'lza', 'm3/mmol N/d', 'quadratic zooplankton loss to dissolved matter', default=0.0666666666_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%lzd, 'lzd', 'm3/mmol N/d', 'quadratic zooplankton loss to detritus', default=0.1333333333_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%iv, 'iv', 'm3/mmol N', 'Ivlev constant for grazing', default=0.24444444_rk)
      call self%get_parameter(self%topt, 'topt', 'deg C', 'optimum temperature for grazing', default=20.0_rk)
      call self%get_parameter(self%lan, 'lan', '1/d', 'maximum nitrification rate', default=0.1_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%oan, 'oan', 'mmol o2/m3', 'half-saturation oxygen concentration for nitrification', default=0.01_rk)
      call self%get_parameter(self%beta_an, 'beta_an', '1/deg C', 'temperature sensitivity of nitrification', default=0.11_rk)
      call self%get_parameter(self%lda, 'lda', '1/d', 'remineralisation rate', default=0.003_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%tda, 'tda', 'deg C', 'half-saturation temperature for remineralisation', default=13.0_rk)
      call self%get_parameter(self%beta_da, 'beta_da', '-', 'temperature sensitivity of remineralisation', default=20.0_rk)
      call self%get_parameter(self%pvel, 'pvel', 'm/d', 'piston velocity', default=5.0_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(self%sr, 'sr', 'mol P/mol N', 'phosphorus to nitrogen ratio', default=0.0625_rk)
      call self%get_parameter(self%s1, 's1', '-', 'reduced nitrate : oxidized detritus for denitrification', default=5.3_rk)
      call self%get_parameter(self%s2, 's2', 'mol O2/mol N', 'oxygen : nitrogen ratio for ammonium-based biomass synthesis and remineralisation', default=6.625_rk)
      call self%get_parameter(self%s3, 's3', 'mol O2/mol N', 'oxygen : nitrogen ratio for nitrate-based biomass synthesis', default=8.625_rk)
      call self%get_parameter(self%s4, 's4', 'mol O2/mol N', 'oxygen : nitrogen ratio for nitrification', default=2.0_rk)
      call self%get_parameter(self%a0, 'a0', 'mmol O2/g', 'mmol O2 per gram', default=31.25_rk)
      call self%get_parameter(self%a1, 'a1', 'g/m3', 'saturation mass concentration of oxygen at 0 deg C', default=14.603_rk)
      call self%get_parameter(self%a2, 'a2', 'g/m3/deg C', 'decrease in saturation mass concentration of oxygen per deg C', default=0.4025_rk)
      if (self%fluff) then
         call self%get_parameter(self%lds, 'lds', 'm/d', 'sedimentation rate of detritus', default=3.5_rk, scale_factor=1.0_rk/secs_pr_day)
         call self%get_parameter(self%lsd, 'lsd', '1/d', 'resuspension rate of fluff', default=25.0_rk, scale_factor=1.0_rk/secs_pr_day)
         call self%get_parameter(self%tau_crit, 'tau_crit', 'critical bottom stress', 'N/m2', default=0.07_rk)
         call self%get_parameter(self%lsa, 'lsa', '1/d', 'fluff mineralisation rate', default=0.001_rk, scale_factor=1.0_rk/secs_pr_day)
         call self%get_parameter(self%bsa, 'bsa', '1/deg C', 'temperature sensitivity of fluff mineralisation', default=0.15_rk)
         call self%get_parameter(self%ph1, 'ph1', '-', 'inhibition of phosphate release during fluff mineralisation', default=0.15_rk)
         call self%get_parameter(self%ph2, 'ph2', 'mmol O2/m3', 'half-saturation oxygen concentration for inhibition of phosphate release', default=0.1_rk)
      end if
      call self%get_parameter(w_p1, 'w_p1', 'vertical velocity of diatoms (< for sinking)', 'm/d', default=-1.0_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(w_p2, 'w_p2', 'vertical velocity of flagellates (< for sinking)', 'm/d', default=-5.0_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(w_p3, 'w_p3', 'vertical velocity of cyanobacteria (< for sinking)', 'm/d', default=-5.0_rk, scale_factor=1.0_rk/secs_pr_day)
      call self%get_parameter(w_de, 'w_de', 'vertical velocity of detritus (< for sinking)', 'm/d', default=-3.0_rk, scale_factor=1.0_rk/secs_pr_day)

      ! Register state variables
      call self%register_state_variable(self%id_p1, 'dia', 'mmol N/m3', 'diatoms',      &
            4.5_rk, minimum=0.0_rk, vertical_movement=w_p1)
      call self%register_state_variable(self%id_p2, 'fla', 'mmol N/m3', 'flagellates',  &
            4.5_rk, minimum=0.0_rk, vertical_movement=w_p2)
      call self%register_state_variable(self%id_p3, 'cya', 'mmol N/m3', 'cyanobacteria',&
            4.5_rk, minimum=0.0_rk, vertical_movement=w_p3)
      call self%register_state_variable(self%id_zo, 'zoo', 'mmol N/m3', 'zooplankton',  &
            4.5_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_de, 'det', 'mmol N/m3', 'detritus',     &
            4.5_rk, minimum=0.0_rk, vertical_movement=w_de)
      call self%register_state_variable(self%id_am, 'amm', 'mmol N/m3', 'ammonium',     &
            4.5_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%register_state_variable(self%id_ni, 'nit', 'mmol N/m3', 'nitrate',      &
            4.5_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%register_state_variable(self%id_po, 'pho', 'mmol P/m3', 'phosphate',    &
            4.5_rk, minimum=0.0_rk, no_river_dilution=.true.)
      call self%register_state_variable(self%id_o2, 'oxy', 'mmol O2/m3', 'oxygen', 4.5_rk)
      if (self%fluff) call self%register_state_variable(self%id_fl, 'flf', 'mmol N/m2', 'fluff', &
            0.0_rk, minimum=0.0_rk)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_dPAR, 'PAR', 'W/m2', 'photosynthetically active radiation')
      call self%register_diagnostic_variable(self%id_GPP,  'GPP', 'mmol/m3', 'gross primary production')
      call self%register_diagnostic_variable(self%id_NCP,  'NCP', 'mmol/m3', 'net community production')
      call self%register_diagnostic_variable(self%id_PPR,  'PPR', 'mmol/m3/d', 'gross primary production rate')
      call self%register_diagnostic_variable(self%id_NPR,  'NPR', 'mmol/m3/d', 'net community production rate')

      ! Register environmental dependencies
      call self%register_dependency(self%id_par,  standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
      call self%register_dependency(self%id_I_0,  standard_variables%surface_downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_wind, standard_variables%wind_speed)
      if (self%fluff) call self%register_dependency(self%id_taub, standard_variables%bottom_stress)

      ! Let phytoplankton (including background concentration) and detritus contribute to light attentuation
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_p1, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_p2, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_p3, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_de, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, (self%p10 + self%p20 + self%p30) * self%kc)
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_gotm_ergom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: p1,p2,p3,zo,am,ni,po,de,o2,par,I_0
      real(rk) :: iopt,ppi,temp,psum,llda,llan,lp,r1,r2,r3
      real(rk) :: wo=30.0_rk,wn=0.1_rk
      real(rk) :: thopnp,thomnp,thomnm,thsum
      real(rk) :: secs_pr_day=86400.0_rk

      ! Enter spatial_loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values
         _GET_(self%id_p1,p1) ! diatoms
         _GET_(self%id_p2,p2) ! flagellates
         _GET_(self%id_p3,p3) ! cyanobacteria
         _GET_(self%id_zo,zo) ! zooplankton
         _GET_(self%id_de,de) ! detritus
         _GET_(self%id_am,am) ! ammonium
         _GET_(self%id_ni,ni) ! nitrate
         _GET_(self%id_po,po) ! phosphate
         _GET_(self%id_o2,o2) ! oxygen

         ! Retrieve current environmental conditions
         _GET_(self%id_par,par)         ! local photosynthetically active radiation
         _GET_SURFACE_(self%id_I_0,I_0) ! surface photosynthetically active radiation
         _GET_(self%id_temp,temp)       ! local water temperature

         ! Light acclimation formulation based on surface light intensity
         iopt = max(0.25*I_0,self%i_min)
         ppi = par/iopt*exp(1.0_rk-par/iopt)

         thopnp=th( o2,wo,0.0_rk,1.0_rk)*yy(wn,ni)
         thomnp=th(-o2,wo,0.0_rk,1.0_rk)*yy(wn,ni)
         thomnm=th(-o2,wo,0.0_rk,1.0_rk)*(1.0_rk-yy(wn,ni))
         thsum=thopnp+thomnp+thomnm
         thopnp=thopnp/thsum
         thomnp=thomnp/thsum
         thomnm=thomnm/thsum

         psum=p1+p2+p3+self%p10+self%p20+self%p30

         llda=self%lda*(1.0_rk+self%beta_da*yy(self%tda,temp))
         llan=th(o2,0.0_rk,0.0_rk,1.0_rk)*o2/(self%oan+o2)*self%lan*exp(self%beta_an*temp)

         lp=self%lpa+self%lpd

         r1=self%r1max*min(yy(self%alpha1,am+ni),yy(self%sr*self%alpha1,po),ppi)*(p1+self%p10)
         r2=self%r2max*(1.0_rk+yy(self%tf,temp))*min(yy(self%alpha2,am+ni),   &
            yy(self%sr*self%alpha2,po),ppi)*(p2+self%p20)
         r3=self%r3max*1.0_rk/(1.0_rk+exp(self%beta_bg*(self%tbg-temp)))       &
            *min(yy(self%sr*self%alpha3,po),ppi)*(p3+self%p30)

         _ADD_SOURCE_(self%id_p1,r1-fpz(self%iv,self%g1max,temp,self%topt,psum)*p1/psum*(zo+self%zo0)-lp*p1)
         _ADD_SOURCE_(self%id_p2,r2-fpz(self%iv,self%g2max,temp,self%topt,psum)*p2/psum*(zo+self%zo0)-lp*p2)
         _ADD_SOURCE_(self%id_p3,r3-fpz(self%iv,self%g3max,temp,self%topt,psum)*p3/psum*(zo+self%zo0)-lp*p3)
         _ADD_SOURCE_(self%id_zo,(fpz(self%iv,self%g1max,temp,self%topt,psum)*p1+fpz(self%iv,self%g2max,temp,self%topt,psum)*p2+fpz(self%iv,self%g3max,temp,self%topt,psum)*p3)*(zo+self%zo0)/psum-self%lza*zo*(zo+self%zo0)-self%lzd*zo*(zo+self%zo0))
         _ADD_SOURCE_(self%id_de,self%lpd*p1+self%lpd*p2+self%lpd*p3+self%lzd*zo*(zo+self%zo0)-llda*de)
         _ADD_SOURCE_(self%id_am,llda*de-llan*am-am/(am+ni)*r1-am/(am+ni)*r2+self%lpa*(p1+p2+p3)+self%lza*zo*(zo+self%zo0))
         _ADD_SOURCE_(self%id_ni,llan*am-ni/(am+ni)*r1-ni/(am+ni)*r2-self%s1*llda*de*thomnp)
         _ADD_SOURCE_(self%id_po,self%sr*(-r1-r2-r3+llda*de+self%lpa*(p1+p2+p3)+self%lza*zo*(zo+self%zo0)))
         _ADD_SOURCE_(self%id_o2,self%s2*(am/(am+ni)*(r1+r2)+r3)-self%s4*(llan*am)+self%s3*(ni/(am+ni)*(r1+r2))-self%s2*(thopnp+thomnm)*llda*de-self%s2*(self%lpa*(p1+p2+p3)+self%lza*zo*(zo+self%zo0)))

         ! Set diagnostic variables
         _SET_DIAGNOSTIC_(self%id_dPAR,par)
         _SET_DIAGNOSTIC_(self%id_GPP ,r1+r2+r3)
         _SET_DIAGNOSTIC_(self%id_NCP ,r1+r2+r3 - self%lpa*(p1+p2+p3))
         _SET_DIAGNOSTIC_(self%id_PPR ,(r1+r2+r3)*secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_NPR ,(r1+r2+r3 - self%lpa*(p1+p2+p3))*secs_pr_day)

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_gotm_ergom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: fl,amb,nib,pob,deb,oxb,taub,temp
      real(rk) :: llds,llsd,llsa,wo=30.0_rk,wn=0.1_rk,dot2=0.2_rk
      real(rk) :: thopnp,thomnp,thomnm,thsum

      if (.not. self%fluff) return

      ! Enter spatial loops over the horizontal domain (if any).
      _BOTTOM_LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_am,amb)
         _GET_(self%id_de,deb)
         _GET_(self%id_ni,nib)
         _GET_(self%id_po,pob)
         _GET_(self%id_o2,oxb)
         _GET_BOTTOM_(self%id_fl,fl)

         _GET_BOTTOM_(self%id_taub,taub)
         _GET_(self%id_temp,temp)

         thopnp=th( oxb,wo,0.0_rk,1.0_rk)*yy(wn,nib)
         thomnp=th(-oxb,wo,0.0_rk,1.0_rk)*yy(wn,nib)
         thomnm=th(-oxb,wo,0.0_rk,1.0_rk)*(1.0_rk-yy(wn,nib))
         thsum=thopnp+thomnp+thomnm
         thopnp=thopnp/thsum
         thomnp=thomnp/thsum
         thomnm=thomnm/thsum

         llsa=self%lsa*exp(self%bsa*temp)*(th(oxb,wo,dot2,1.0_rk))

         ! Sedimentation dependent on bottom stress
         if (self%tau_crit > taub) then
            llds=self%lds*(self%tau_crit-taub)/self%tau_crit
         else
            llds=0.
         end if

         ! Resuspension dependent on bottom stress
         if (self%tau_crit < taub) then
            llsd=self%lsd*(taub-self%tau_crit)/self%tau_crit
         else
            llsd=0.0_rk
         end if

         _ADD_BOTTOM_SOURCE_(self%id_fl,llds*deb-llsd*fl-llsa*fl-th(oxb,wo,0.0_rk,1.0_rk)*llsa*fl)
         _ADD_BOTTOM_FLUX_(self%id_de,-llds*deb+llsd*fl)
         _ADD_BOTTOM_FLUX_(self%id_am,llsa*fl)
         _ADD_BOTTOM_FLUX_(self%id_ni,-self%s1*thomnp*llsa*fl)
         _ADD_BOTTOM_FLUX_(self%id_po,self%sr*(1.0_rk-self%ph1*th(oxb,wo,0.0_rk,1.0_rk)*yy(self%ph2,oxb))*llsa*fl)
         _ADD_BOTTOM_FLUX_(self%id_o2,-(self%s4+self%s2*(thopnp+thomnm))*llsa*fl)

      ! Leave spatial loops over the horizontal domain (if any).
      _BOTTOM_LOOP_END_
   end subroutine do_bottom

   ! Weiss formula for the saturation oxygen (osat) \cite{Weiss1970}:
   !
   ! \begin{equation}\label{osat_weiss}
   ! O_{sat}= \exp\left[a_1 +a_2\frac{100}{T}+a_3\ln\left(\frac{T}{100}\right)
   ! +a_4\frac{T}{100}+S \left\{b_1+b_2\frac{T}{100}
   ! +b_3\left(\frac{T}{100}\right)^2   \right\}\right],
   ! \end{equation}
   !
   ! where $T$ is the temperature in Kelvin and the empirical constants are
   ! given in table \ref{table_weiss}.
   !
   ! \begin{table}[h]
   ! \begin{center}
   ! \begin{tabular}{|l|l|l|l|l|l|l|}
   ! \hline
   ! $a_1$ & $a_2$ & $a_3$ & $a_4$ & $b_1$ & $b_2$ & $b_3$ \\ \hline
   ! -173.4292 &
   ! 249.6339 &
   ! 143.3483 &
   ! -21.8492 &
   ! -0.033096 &
   ! 0.014259 &
   ! -0.001700 \\ \hline
   ! \end{tabular}
   ! \caption{Constants for the oxygen saturation formula by \cite{Weiss1970},
   ! see equation (\ref{osat_weiss}).}
   ! \label{table_weiss}
   ! \end{center}
   ! \end{table}
   elemental real(rk) function osat_weiss(t,s)
      real(rk), intent(in) :: t, s

      real(rk)            :: tk
      real(rk), parameter :: aa1=-173.4292_rk
      real(rk), parameter :: aa2=249.6339_rk
      real(rk), parameter :: a3=143.3483_rk
      real(rk), parameter :: a4=-21.8492_rk
      real(rk), parameter :: b1=-0.033096_rk
      real(rk), parameter :: b2=0.014259_rk
      real(rk), parameter :: b3=-0.001700_rk
      real(rk), parameter :: kelvin=273.16_rk
      real(rk), parameter :: mol_per_liter=44.661_rk

      tk=(t+kelvin)*0.01_rk
      osat_weiss=exp(aa1+aa2/tk+a3*log(tk)+a4*tk    &
                 +s*(b1+(b2+b3*tk)*tk))*mol_per_liter
   end function osat_weiss

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_gotm_ergom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk)            :: temp,wnd,salt,o2,ni,am,po
      real(rk), parameter :: secs_pr_day=86400.0_rk
      real(rk)            :: p_vel,sc,flo2
      integer, parameter  :: newflux=1

      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_temp,temp)
         _GET_(self%id_salt,salt)
         _GET_SURFACE_(self%id_wind,wnd)

         _GET_(self%id_o2,o2)
         _GET_(self%id_ni,ni)
         _GET_(self%id_am,am)
         _GET_(self%id_po,po)

         ! Calculation of the surface oxygen flux
         if (newflux .eq. 1) then
            sc=1450.+(1.1*temp-71.0_rk)*temp
            if (wnd .gt. 13.0_rk) then
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
            p_vel = p_vel/secs_pr_day
            flo2 =p_vel*(osat_weiss(temp,salt)-o2)
         else
            flo2 = self%pvel*(self%a0*(self%a1-self%a2*temp)-o2)
         end if
         _ADD_SURFACE_FLUX_(self%id_o2,flo2)

         _ADD_SURFACE_FLUX_(self%id_ni,self%sfl_ni)
         _ADD_SURFACE_FLUX_(self%id_am,self%sfl_am)
         _ADD_SURFACE_FLUX_(self%id_po,self%sfl_po)

      _SURFACE_LOOP_END_
   end subroutine do_surface

   ! Step function
   !
   ! Instead of the
   ! heavyside switches used by \cite{Neumannetal2002}, we apply here a smoothed
   ! {\it tangens hyperbolicus} transition with prescribed width $x_w$:
   ! \begin{equation}\label{theta}
   ! \theta (x,x_w,y_{\min},y_{\max})= y_{\min}+(y_{\max}-y_{\min})
   ! \frac12\left(1-\tanh \left(\frac{x}{x_w}   \right)      \right).
   ! \end{equation}
   elemental real(rk) function th(x,w,min,max)
      real(rk), intent(in) :: x,w,min,max

      if (w .gt. 1.e-10_rk) then
         th=min+(max-min)*0.5_rk*(1.0_rk+tanh(x/w))
      else
         if(x .gt. 0.0_rk) then
            th=1.0_rk
         else
            th=0.0_rk
         end if
      end if
   end function th

   ! Functional response for nutrient uptake
   !
   ! This is a squared Michaelis-Menten type of limiter:
   ! \begin{equation}\label{Y}
   ! Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
   ! \end{equation}
   elemental real(rk) function yy(a,x)
      real(rk), intent(in) :: a,x

      yy=x**2/(a**2+x**2)
   end function yy

   ! Ivlev formulation for zooplankton grazing on phytoplankton
   !
   ! The Ivlev formulation for zooplankton grazing on the three phytoplankton
   ! classes $c_1$, $c_2$, and $c_3$ is given here as a function:
   ! \begin{equation}\label{neu_di4}
   ! d_{i,4}=g_i^{\max}\left(1+\frac{T^2}{T_{opt}^2}\exp
   ! \left(1-\frac{2T}{T_{opt}} \right)\right)
   ! \left( 1-\exp\left(-I_v^2 \left( \sum_{
   ! j=1}^3c_j \right)^2\right)  \right)
   ! \frac{c_i}{\sum_{j=1}^3c_j}\left( c_4+c_4^{\min} \right)
   ! \end{equation}
   pure real(rk) function fpz(iv,g,t,topt,psum)
      real(rk), intent(in) :: iv,g,t,topt,psum

      fpz = g * (1.0_rk+t**2/topt**2*exp(1.0_rk-2.0_rk*t/topt)) * (1.0_rk-exp(-iv**2*psum**2))
   end function fpz

  end module gotm_ergom

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

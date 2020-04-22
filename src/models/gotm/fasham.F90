#include "fabm_driver.h"

! Fasham et al. (1990, DSR I) biogeochemical model,
! with slight modifications by Kuehn \& Radach (1997, JRM),
! and minimum concentration for phytoplankton, zooplankton and bacteria
! as decribed in Burchard et al. (2006, JMS)
! taken from GOTM and adapted for FABM by Jorn Bruggeman
!
!  The model developed by \cite{Fashametal1990} 
!  uses nitrogen as 'currency' according to the evidence that in
!  most cases nitrogen is the limiting macronutrient. It consists of
!  seven state variables: phytoplankton, zooplankton, bacteria,
!  particulate organic matter (detritus), dissolved organic matter
!  and the nutrients nitrate and ammonium.
!  The structure of the \cite{Fashametal1990} biogeochemical model
!  is given in figure \ref{fig_fasham}.
! \begin{figure}
! \begin{center}
! \scalebox{0.5}{\includegraphics{figures/fasham_structure.eps}}
! \caption{Structure of the \cite{Fashametal1990} model with bacteria (bac),
! phytoplankton (phy), detritus (det), zooplankton (zoo), labile dissolved
! organic nitrogen (don), ammonium (amm) and nitrate (nit) as the seven
! state variables.
! The concentrations are in mmol N\,m$^{-3}$,
! all fluxes (green arrows) are conservative.
! }\label{fig_fasham}
! \end{center}
! \end{figure}
!  A detailed mathematical description of all
!  processes is given in section \ref{sec:bio-fasham-rhs}.
!  The version of the \cite{Fashametal1990} model which is implemented includes
!  slight modifications by \cite{KuehnRadach1997} and has been 
!  included into GOTM by \cite{Burchardetal05}.
!
! The \cite{Fashametal1990} model consisting of the $I=7$
! state variables phytoplankton, bacteria, detritus, zooplankton, 
! nitrate, ammonium and dissolved organic nitrogen is described here
! in detail.
! 
! Phytoplankton mortality and zooplankton grazing loss of phytoplankton:
! \begin{equation}\label{d13}
! d_{1,3} = \mu_1 \frac{c_1+c_{1}^{\min}}{K_5+c_1+c_{1}^{\min}}c_1+
! (1-\beta)\frac{g\rho_1 c_1^2}{K_3 \sum_{j=1}^3 \rho_jc_j
! + \sum_{j=1}^3 \rho_jc_j^2} (c_4+c_{4}^{\min}).
! \end{equation}
! Phytoplankton loss to LDON (labile dissolved organic nitrogen):
! \begin{equation}\label{d17}
! d_{1,7} = \gamma
! F(I_{PAR})\frac{\frac{c_5}{K_1}
! +\frac{c_6}{K_2}}{1+\frac{c_5}{K_1}+\frac{c_6}{K_2}}c_1,
! \end{equation}
! with
! \begin{equation}\label{FI}
!  F(I_{PAR}) = \frac{V_p\alpha I_{PAR}(z)}{\left(V_p^2+\alpha^2(I_{PAR}(z))^2 
! \right)^{1/2}}.
! \end{equation}
! With $I_{PAR}$ from (\ref{light}). 
! 
! Zooplankton grazing loss:
! \begin{equation}\label{di3}
! d_{2,3} = (1-\beta)\frac{g\rho_2 c_2^2}{K_3 \sum_{j=1}^3 \rho_jc_j 
! + \sum_{j=1}^3 \rho_jc_j^2} (c_4+c_{4}^{\min}).
! \end{equation}
! Zooplankton grazing:
! \begin{equation}\label{di4}
! d_{i,4} = \beta\frac{g\rho_i c_i^2}{K_3 \sum_{j=1}^3 \rho_jc_j 
! + \sum_{j=1}^3 \rho_jc_j^2} (c_4+c_{4}^{\min}), \quad i=1,\dots,3.
! \end{equation}
! Bacteria excretion rate:
! \begin{equation}\label{d26}
! d_{2,6} = \mu_3 c_2.
! \end{equation}
! Detritus breakdown rate:
! \begin{equation}\label{d37}
! d_{3,7} = \mu_4 c_3.
! \end{equation}
! Zooplankton losses to detritus, ammonium and LDON:
! \begin{equation}\label{d43}
! d_{4,3} = (1-\epsilon-\delta)\mu_2 
! \frac{c_4+c_{4}^{\min}}{K_6+c_4+c_{4}^{\min}}c_4.
! \end{equation}
! \begin{equation}\label{d46}
! d_{4,6} = \epsilon\mu_2 \frac{c_4+c_{4}^{\min}}{K_6+c_4+c_{4}^{\min}}c_4.
! \end{equation}
! \begin{equation}\label{d47}
! d_{4,7} = \delta\mu_2 \frac{c_4+c_{4}^{\min}}{K_6+c_4+c_{4}^{\min}}c_4.
! \end{equation}
! Nitrate uptake by phytoplankton:
! \begin{equation}\label{d51}
! d_{5,1} = F(I_{PAR})\frac{\frac{c_5}{K_1}}{1+\frac{c_5}{K_1}
! +\frac{c_6}{K_2}}(c_1+c_{1}^{\min}).
! \end{equation}
! Ammonium uptake by phytoplankton:
! \begin{equation}\label{d61}
! d_{6,1} = F(I_{PAR})\frac{\frac{c_6}{K_2}}{1+\frac{c_5}{K_1}
! +\frac{c_6}{K_2}}(c_1+c_{1}^{\min}).
! \end{equation}
! Ammonium uptake by bacteria:
! \begin{equation}\label{d62}
! d_{6,2} = V_b \frac{\min(c_6,\eta c_7)}{K_4+\min(c_6,\eta c_7)+c_7} 
! (c_2+c_{2}^{\min}).
! \end{equation}
! LDON uptake by bacteria:
! \begin{equation}\label{d72}
! d_{7,2} = V_b \frac{c_7}{K_4+\min(c_6,\eta c_7)+c_7} (c_2+c_{2}^{\min}).
! \end{equation}

module gotm_fasham

   use fabm_types
   
   implicit none

   private

   type, extends(type_base_model), public :: type_gotm_fasham
      ! Variable identifiers
      type (type_state_variable_id)      :: id_p, id_z, id_b, id_d, id_n, id_a, id_l
      type (type_diagnostic_variable_id) :: id_pp
      type (type_dependency_id)          :: id_par
      
      ! Model parameters
      real(rk) :: p0
      real(rk) :: z0
      real(rk) :: b0
      real(rk) :: vp
      real(rk) :: alpha
      real(rk) :: k1
      real(rk) :: k2
      real(rk) :: mu1
      real(rk) :: k5
      real(rk) :: gamma
      real(rk) :: gmax
      real(rk) :: k3
      real(rk) :: beta
      real(rk) :: mu2
      real(rk) :: k6
      real(rk) :: delta
      real(rk) :: epsi
      real(rk) :: r1
      real(rk) :: r2
      real(rk) :: r3
      real(rk) :: vb
      real(rk) :: k4
      real(rk) :: mu3
      real(rk) :: eta
      real(rk) :: mu4
      real(rk) :: kc
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
   end type

contains

   subroutine initialize(self, configunit)
      class (type_gotm_fasham), intent(inout),target :: self
      integer,                  intent(in)           :: configunit

      real(rk), parameter :: days_per_sec = 1.0_rk/86400.0_rk
      real(rk)            :: w_p, w_d

      ! Store parameter values in our own derived type
      call self%get_parameter(self%p0,   'p0',   'mmol m-3',  'minimum phytoplankton concentration',default=0.0_rk)
      call self%get_parameter(self%vp,   'vp',   'd-1',       'maximum phytoplankton uptake rate',default=1.5_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%alpha,'alpha','m2 W-1 d-1','initial slope of photosynthesis-irradiance curve',default=0.065_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%k1,   'k1',   'mmol m-3',  'half saturation constant for nitrate uptake',default=0.2_rk)
      call self%get_parameter(self%k2,   'k2',   'mmol m-3',  'half saturation constant for ammonium uptake',default=0.8_rk)
      call self%get_parameter(self%mu1,  'mu1',  'd-1',       'maximum phytoplankton mortality rate',default=0.05_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%k5,   'k5',   'mmol m-3',  'half saturation constant for phytoplankton mortality',default=0.2_rk)
      call self%get_parameter(self%gamma,'gamma','-',         'fraction of primary production that is exudated',default=0.05_rk)
      call self%get_parameter(w_p,       'w_p',  'm d-1',     'phytoplankton settling velocity (negative for sinking)',default=-1.0_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%kc,   'kc',   'm2 mmol-1', 'specific light attenuation of phytoplankton',default=0.03_rk)
      call self%get_parameter(self%z0,   'z0',   'mmol m-3',  'minimum zooplankton concentration',default=0.0_rk)
      call self%get_parameter(self%gmax, 'gmax', 'd-1',       'maximum ingestion rate',default=1.0_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%k3,   'k3',   'mmol m-3',  'half saturation constant for zooplankton ingestion',default=1.0_rk)
      call self%get_parameter(self%beta, 'beta', '-',         'grazing efficiency',default=0.625_rk)
      call self%get_parameter(self%mu2,  'mu2',  'd-1',       'maximum zooplankton loss rate',default=0.3_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%k6,   'k6',   'mmol m-3',  'half saturation constant for zooplankton loss',default=0.2_rk)
      call self%get_parameter(self%delta,'delta','-',         'fractional zooplankton loss to LDON',default=0.1_rk)
      call self%get_parameter(self%epsi, 'epsi', '-',         'fractional zooplankton loss to ammonium',default=0.7_rk)
      call self%get_parameter(self%r1,   'r1',   '-',         'zooplankton preference for phytoplankton',default=0.55_rk)
      call self%get_parameter(self%r2,   'r2',   '-',         'zooplankton preference for bacteria',default=0.4_rk)
      call self%get_parameter(self%r3,   'r3',   '-',         'zooplankton preference for detritus',default=0.05_rk)
      call self%get_parameter(self%b0,   'b0',   'mmol m-3',  'minimum bacteria concentration',default=0.0_rk)
      call self%get_parameter(self%vb,   'vb',   'd-1',       'maximum bacterial uptake rate',default=1.2_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%k4,   'k4',   'mmol m-3',  'half saturation constant for bacterial uptake',default=0.5_rk)
      call self%get_parameter(self%mu3,  'mu3',  'd-1',       'bacterial excretion rate',default=0.15_rk, scale_factor=days_per_sec)
      call self%get_parameter(self%eta,  'eta',  '-',         'bacterial ammonium:LDON uptake ratio',default=0.0_rk)
      call self%get_parameter(self%mu4,  'mu4',  'd-1',       'detritus breakdown rate',default=0.02_rk, scale_factor=days_per_sec)
      call self%get_parameter(w_d,       'w_d',  'm d-1',     'detritus settling velocity (negative for sinking)',default=-2.0_rk, scale_factor=days_per_sec)

      ! Register state variables
      call self%register_state_variable(self%id_p,'phy','mmol m-3','phytoplankton',     &
                                       0.056666666_rk,minimum=0.0_rk,vertical_movement=w_p)
      call self%register_state_variable(self%id_z,'zoo','mmol m-3','zooplankton',     &
                                       0.05_rk,minimum=0.0_rk)
      call self%register_state_variable(self%id_b,'bac','mmol m-3','bacteria',     &
                                       0.001_rk,minimum=0.0_rk)
      call self%register_state_variable(self%id_d,'det','mmol m-3','detritus',     &
                                       0.416666666_rk,minimum=0.0_rk,vertical_movement=w_d)
      call self%register_state_variable(self%id_n,'nit','mmol m-3','nitrate',     &
                                       8.3_rk,minimum=0.0_rk,no_river_dilution=.true.)
      call self%register_state_variable(self%id_a,'amm','mmol m-3','ammonium',     &
                                       0.22_rk,minimum=0.0_rk,no_river_dilution=.true.)
      call self%register_state_variable(self%id_l,'ldn','mmol m-3','labile dissolved organic nitrogen',     &
                                       0.14_rk,minimum=0.0_rk,no_river_dilution=.true.)

      ! Register the contribution of all state variables to total nitrogen
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_p)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_z)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_b)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_d)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_a)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_l)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_pp,'pp','d-1','specific primary production rate')

      ! Register environmental dependencies
      call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)

      ! Let phytoplankton (incldung background concentration) contribute to light attenuation
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_p, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%p0 * self%kc)
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_gotm_fasham), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: p,z,b,d,n,a,l,par
      real(rk)            :: d_p,d_z,d_b,d_d,d_n,d_a,d_l
      real(rk)            :: ff,fac,min67
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_p,p) ! phytoplankton
         _GET_(self%id_z,z) ! zooplankton
         _GET_(self%id_b,b) ! bacteria
         _GET_(self%id_d,d) ! detritus
         _GET_(self%id_n,n) ! nitrate
         _GET_(self%id_a,a) ! ammonia
         _GET_(self%id_l,l) ! labile dissolved organic nitrogen
   
         ! Retrieve local environment: photosynthetically active radiation.
         _GET_(self%id_par,par)

         ! Calculate intermediate quantities:
         !   ff    = phytoplankton light limitation
         !   fac   = zooplankton grazing denominator multiplied by zooplankton
         !   min67 = total bacterial nitrogenous substrate
         ff = self%vp*self%alpha*par/sqrt(self%vp**2+self%alpha**2*par**2) 
         if (p.eq.0.0_rk .and. b.eq.0.0_rk .and. d.eq.0.0_rk) then
            fac = 0.0_rk
         else
            fac = (z+self%z0)/(self%k3*(self%r1*p+self%r2*b+self%r3*d)+  &
                               self%r1*p**2+self%r2*b**2+self%r3*d**2)
         end if
         min67 = min(a,self%eta*l)

         ! Calculate temporal derivatives according to Kuehn & Radach (1997, Journal of Marine Research)
         d_p = ff*(n/self%k1+a/self%k2)/(1.0_rk+n/self%k1+a/self%k2)*(p*(1.0_rk-self%gamma)+self%p0) &
               - self%gmax*self%r1*p**2*fac                                                        &
               - self%mu1*(p+self%p0)/(self%k5+p+self%p0)*p
         d_z =   self%beta*self%gmax*self%r1*p**2*fac &
               + self%beta*self%gmax*self%r2*b**2*fac &
               + self%beta*self%gmax*self%r3*d**2*fac &
               - self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z
         d_b =   self%vb*l    /(self%k4+min67+l)*(b+self%b0) &
               + self%vb*min67/(self%k4+min67+l)*(b+self%b0) &
               - self%gmax*self%r2*b**2*fac                  &
               - self%mu3*b
         d_d =   (1.0_rk-self%beta)*self%gmax*self%r1*p**2*fac &
               + (1.0_rk-self%beta)*self%gmax*self%r2*b**2*fac &
               - self%beta        *self%gmax*self%r3*d**2*fac &
               - self%mu4*d                                   &
               + self%mu1*(p+self%p0)/(self%k5+p+self%p0)*p   &
               + (1.0_rk-self%epsi-self%delta)*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z
         d_n = - ff*n/self%k1/(1.0_rk+n/self%k1+a/self%k2)*(p+self%p0)
         d_a = - ff*a/self%k2/(1.0_rk+n/self%k1+a/self%k2)*(p+self%p0) &
               - self%vb*min67/(self%k4+min67+l)*(b+self%b0)          &
               + self%mu3*b                                           &
               + self%epsi*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z
         d_l =   self%gamma*ff*(n/self%k1+a/self%k2)/(1.0_rk+n/self%k1+a/self%k2)*p &
               + self%mu4*d                                                        &
               + self%delta*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z             &
               - self%vb*l/(self%k4+min67+l)*(b+self%b0)

         ! Provide temporal derivatives to FABM.
         _ADD_SOURCE_(self%id_p,d_p)
         _ADD_SOURCE_(self%id_z,d_z)
         _ADD_SOURCE_(self%id_d,d_d)
         _ADD_SOURCE_(self%id_b,d_b)
         _ADD_SOURCE_(self%id_n,d_n)
         _ADD_SOURCE_(self%id_a,d_a)
         _ADD_SOURCE_(self%id_l,d_l)
   
         ! Provide diagnostic variables to FABM.
         _SET_DIAGNOSTIC_(self%id_pp,secs_pr_day*ff*(n/self%k1+a/self%k2)/(1.0_rk+n/self%k1+a/self%k2))
   
      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

   subroutine do_ppdd(self, _ARGUMENTS_DO_PPDD_)
      class (type_gotm_fasham), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_PPDD_

      real(rk)            :: p,z,b,d,n,a,l,par
      real(rk)            :: ff,fac,min67
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_p,p) ! phytoplankton
         _GET_(self%id_z,z) ! zooplankton
         _GET_(self%id_b,b) ! bacteria
         _GET_(self%id_d,d) ! detritus
         _GET_(self%id_n,n) ! nitrate
         _GET_(self%id_a,a) ! ammonia
         _GET_(self%id_l,l) ! labile dissolved organic nitrogen
   
         ! Retrieve local photosynthetically active radiation.
         _GET_(self%id_par,par)

         ! Calculate intermediate quantities that will be used multiple times.
         ff= self%vp*self%alpha*par/sqrt(self%vp**2+self%alpha**2*par**2)
         if (p.eq.0.0_rk .and. b.eq.0.0_rk .and. d.eq.0.0_rk) then
            fac = 0.0_rk
         else
           fac=(z+self%z0)/(self%k3*(self%r1*p+self%r2*b+self%r3*d)+  &
                            self%r1*p**2+self%r2*b**2+self%r3*d**2)
         end if
         min67=min(a,self%eta*l)

         ! Set production & destruction terms. Note that SET_DD_SYM is used to assign a
         ! destruction term, but it will automatically set the corresponding production
         ! term as well, i.e., pp(j,i) = dd(i,j).
         _SET_DD_SYM_(self%id_p,self%id_d,self%mu1*(p+self%p0)/(self%k5+p+self%p0)*p+(1.0_rk-self%beta)*self%gmax*self%r1*p**2*fac)
         _SET_DD_SYM_(self%id_p,self%id_l,self%gamma*ff*(n/self%k1+a/self%k2)/(1.0_rk+n/self%k1+a/self%k2)*p)
         _SET_DD_SYM_(self%id_b,self%id_d,(1.0_rk-self%beta)*self%gmax*self%r2*b**2*fac)
         _SET_DD_SYM_(self%id_p,self%id_z,self%beta*self%gmax*self%r1*p**2*fac)
         _SET_DD_SYM_(self%id_b,self%id_z,self%beta*self%gmax*self%r2*b**2*fac)
         _SET_DD_SYM_(self%id_d,self%id_z,self%beta*self%gmax*self%r3*d**2*fac)
         _SET_DD_SYM_(self%id_b,self%id_a,self%mu3*b)
         _SET_DD_SYM_(self%id_d,self%id_l,self%mu4*d)
         _SET_DD_SYM_(self%id_z,self%id_d,(1.0_rk-self%epsi-self%delta)*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z)
         _SET_DD_SYM_(self%id_z,self%id_a,self%epsi*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z)
         _SET_DD_SYM_(self%id_z,self%id_l,self%delta*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z)
         _SET_DD_SYM_(self%id_n,self%id_p,ff*n/self%k1/(1.0_rk+n/self%k1+a/self%k2)*(p+self%p0))
         _SET_DD_SYM_(self%id_a,self%id_p,ff*a/self%k2/(1.0_rk+n/self%k1+a/self%k2)*(p+self%p0))
         _SET_DD_SYM_(self%id_a,self%id_b,self%vb*min67/(self%k4+min67+l)*(b+self%b0))
         _SET_DD_SYM_(self%id_l,self%id_b,self%vb*l/(self%k4+min67+l)*(b+self%b0))
   
         ! Provide diagnostic variables to FABM.
         _SET_DIAGNOSTIC_(self%id_pp,secs_pr_day*ff*(n/self%k1+a/self%k2)/(1.0_rk+n/self%k1+a/self%k2))

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do_ppdd

end module gotm_fasham

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

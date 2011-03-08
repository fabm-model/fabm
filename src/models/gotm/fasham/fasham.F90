!$Id: fasham.F90 119 2010-12-27 14:23:18Z jornbr $
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_fasham --- Fasham et al. (1990, DSR I) biogeochemical model,
! with slight modifications by Kuehn & Radach (1997, JRM),
! taken from GOTM and adapted for FABM by Jorn Bruggeman
!
! !INTERFACE:
   module fabm_fasham
!
! !DESCRIPTION:
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
! !USES:
   use fabm_types
   use fabm_driver

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_fasham, fasham_init, fasham_do, fasham_do_ppdd, &
          fasham_get_light_extinction, fasham_get_conserved_quantities
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!
! !PUBLIC DERIVED TYPES:
   type type_fasham
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_      :: id_p,id_z,id_b,id_d,id_n,id_a,id_l
      _TYPE_DEPENDENCY_ID_          :: id_par
      _TYPE_CONSERVED_QUANTITY_ID_  :: id_totN
      
!     Model parameters
      REALTYPE :: p0
      REALTYPE :: z0
      REALTYPE :: b0
      REALTYPE :: vp
      REALTYPE :: alpha
      REALTYPE :: k1
      REALTYPE :: k2
      REALTYPE :: mu1
      REALTYPE :: k5
      REALTYPE :: gamma
      REALTYPE :: gmax
      REALTYPE :: k3
      REALTYPE :: beta
      REALTYPE :: mu2
      REALTYPE :: k6
      REALTYPE :: delta
      REALTYPE :: epsi
      REALTYPE :: r1
      REALTYPE :: r2
      REALTYPE :: r3
      REALTYPE :: vb
      REALTYPE :: k4
      REALTYPE :: mu3
      REALTYPE :: eta
      REALTYPE :: mu4
      REALTYPE :: kc
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the Fasham model
!
! !INTERFACE:
   subroutine fasham_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the fasham namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_fasham),      intent(out)   :: self
   type (type_model_info),intent(inout) :: modelinfo
   integer,               intent(in)    :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                  ::  p_initial= 0.056666666
   REALTYPE                  ::  z_initial= 0.05
   REALTYPE                  ::  b_initial= 0.001
   REALTYPE                  ::  d_initial= 0.416666666
   REALTYPE                  ::  n_initial= 8.3
   REALTYPE                  ::  a_initial= 0.22
   REALTYPE                  ::  l_initial= 0.14
   REALTYPE                  ::  p0       = 0.0
   REALTYPE                  ::  z0       = 0.0
   REALTYPE                  ::  b0       = 0.0
   REALTYPE                  ::  kc       = 0.03
   REALTYPE                  ::  vp       = 1.5
   REALTYPE                  ::  alpha    = 0.065
   REALTYPE                  ::  k1       = 0.2
   REALTYPE                  ::  k2       = 0.8
   REALTYPE                  ::  mu1      = 0.05
   REALTYPE                  ::  k5       = 0.2
   REALTYPE                  ::  gamma    = 0.05
   REALTYPE                  ::  w_p      = -1.0
   REALTYPE                  ::  gmax     = 1.0
   REALTYPE                  ::  k3       = 1.0
   REALTYPE                  ::  beta     = 0.625
   REALTYPE                  ::  mu2      = 0.3
   REALTYPE                  ::  k6       = 0.2
   REALTYPE                  ::  delta    = 0.1
   REALTYPE                  ::  epsi     = 0.70
   REALTYPE                  ::  r1       = 0.55
   REALTYPE                  ::  r2       = 0.4
   REALTYPE                  ::  r3       = 0.05
   REALTYPE                  ::  vb       = 1.2
   REALTYPE                  ::  k4       = 0.5
   REALTYPE                  ::  mu3      = 0.15
   REALTYPE                  ::  eta      = 0.0
   REALTYPE                  ::  mu4      = 0.02
   REALTYPE                  ::  w_d      = -2.0

   REALTYPE, parameter :: secs_pr_day = 86400.
   namelist /fasham/ p_initial,z_initial,b_initial,d_initial,n_initial,&
                     a_initial,l_initial,p0,z0,b0,vp,alpha,k1,k2,mu1,k5,&
                     gamma,w_p,gmax,k3,beta,mu2,k6,delta,epsi,r1,r2,r3, &
                     vb,k4,mu3,eta,mu4,w_d,kc
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=fasham,err=99,end=100)

   ! All rates must be provided in values per day,
   ! and are converted here to values per second.
   vp   = vp   /secs_pr_day
   vb   = vb   /secs_pr_day
   mu1  = mu1  /secs_pr_day
   mu2  = mu2  /secs_pr_day
   mu3  = mu3  /secs_pr_day
   mu4  = mu4  /secs_pr_day
   gmax = gmax /secs_pr_day
   w_p  = w_p  /secs_pr_day
   w_d  = w_d  /secs_pr_day
   alpha= alpha/secs_pr_day

   ! Store parameter values in our own derived type
   self%p0    = p0
   self%z0    = z0
   self%b0    = b0
   self%vp    = vp
   self%alpha = alpha
   self%k1    = k1
   self%k2    = k2
   self%mu1   = mu1
   self%k5    = k5
   self%gamma = gamma
   self%gmax  = gmax
   self%k3    = k3
   self%beta  = beta
   self%mu2   = mu2
   self%k6    = k6
   self%delta = delta
   self%epsi  = epsi
   self%r1    = r1
   self%r2    = r2
   self%r3    = r3
   self%vb    = vb
   self%k4    = k4
   self%mu3   = mu3
   self%eta   = eta
   self%mu4   = mu4
   self%kc    = kc
   
   ! Register state variables
   self%id_p = register_state_variable(modelinfo,'phy','mmol/m**3','phytoplankton',     &
                                    p_initial,minimum=_ZERO_,vertical_movement=w_p)
   self%id_z = register_state_variable(modelinfo,'zoo','mmol/m**3','zooplankton',     &
                                    z_initial,minimum=_ZERO_)
   self%id_b = register_state_variable(modelinfo,'bac','mmol/m**3','bacteria',     &
                                    b_initial,minimum=_ZERO_)
   self%id_d = register_state_variable(modelinfo,'det','mmol/m**3','detritus',     &
                                    d_initial,minimum=_ZERO_,vertical_movement=w_d)
   self%id_n = register_state_variable(modelinfo,'nit','mmol/m**3','nitrate',     &
                                    n_initial,minimum=_ZERO_,no_river_dilution=.true.)
   self%id_a = register_state_variable(modelinfo,'amm','mmol/m**3','ammonium',     &
                                    a_initial,minimum=_ZERO_,no_river_dilution=.true.)
   self%id_l = register_state_variable(modelinfo,'ldn','mmol/m**3','labile dissolved organic nitrogen',     &
                                    l_initial,minimum=_ZERO_,no_river_dilution=.true.)

   ! Register environmental dependencies
   self%id_par = register_dependency(modelinfo, varname_par)

   ! Register conserved quantities
   self%id_totN = register_conserved_quantity(modelinfo,'N','mmol/m**3','nitrogen')

   return

99 call fatal_error('fasham_init','Error reading namelist fasham')
100 call fatal_error('fasham_init','Namelist fasham was not found')
   
   end subroutine fasham_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Fasham model
!
! !INTERFACE:
   subroutine fasham_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!
! !DESCRIPTION:
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
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_fasham),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                   :: p,z,b,d,n,a,l,par
   REALTYPE                   :: ff,fac,min67
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_p,p) ! phytoplankton
   _GET_STATE_(self%id_z,z) ! zooplankton
   _GET_STATE_(self%id_b,b) ! bacteria
   _GET_STATE_(self%id_d,d) ! detritus
   _GET_STATE_(self%id_n,n) ! nitrate
   _GET_STATE_(self%id_a,a) ! ammonia
   _GET_STATE_(self%id_l,l) ! labile dissolved organic nitrogen
   
   ! Retrieve local photosynthetically active radiation.
   _GET_DEPENDENCY_(self%id_par,par)

   ! Calculate intermediate quantities that will be used multiple times.
   ff= self%vp*self%alpha*par/sqrt(self%vp**2+self%alpha**2*par**2)
   if (p.eq._ZERO_ .and. b.eq._ZERO_ .and. d.eq._ZERO_) then
      fac = _ZERO_
   else
     fac=(z+self%z0)/(self%k3*(self%r1*p+self%r2*b+self%r3*d)+  &
                      self%r1*p**2+self%r2*b**2+self%r3*d**2)
   end if
   min67=min(a,self%eta*l)

   ! Set production & destruction terms. Note that SET_DD_SYM is used to assign a
   ! destruction term, but it will automatically set the corresponding production
   ! term as well, i.e., pp(j,i) = dd(i,j).
   _SET_DD_SYM_(self%id_p,self%id_d,self%mu1*(p+self%p0)/(self%k5+p+self%p0)*p+(_ONE_-self%beta)*self%gmax*self%r1*p**2*fac)
   _SET_DD_SYM_(self%id_p,self%id_l,self%gamma*ff*(n/self%k1+a/self%k2)/(_ONE_+n/self%k1+a/self%k2)*p)
   _SET_DD_SYM_(self%id_b,self%id_d,(_ONE_-self%beta)*self%gmax*self%r2*b**2*fac)
   _SET_DD_SYM_(self%id_p,self%id_z,self%beta*self%gmax*self%r1*p**2*fac)
   _SET_DD_SYM_(self%id_b,self%id_z,self%beta*self%gmax*self%r2*b**2*fac)
   _SET_DD_SYM_(self%id_d,self%id_z,self%beta*self%gmax*self%r3*d**2*fac)
   _SET_DD_SYM_(self%id_b,self%id_a,self%mu3*b)
   _SET_DD_SYM_(self%id_d,self%id_l,self%mu4*d)
   _SET_DD_SYM_(self%id_z,self%id_d,(_ONE_-self%epsi-self%delta)*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z)
   _SET_DD_SYM_(self%id_z,self%id_a,self%epsi*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z)
   _SET_DD_SYM_(self%id_z,self%id_l,self%delta*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z)
   _SET_DD_SYM_(self%id_n,self%id_p,ff*n/self%k1/(_ONE_+n/self%k1+a/self%k2)*(p+self%p0))
   _SET_DD_SYM_(self%id_a,self%id_p,ff*a/self%k2/(_ONE_+n/self%k1+a/self%k2)*(p+self%p0))
   _SET_DD_SYM_(self%id_a,self%id_b,self%vb*min67/(self%k4+min67+l)*(b+self%b0))
   _SET_DD_SYM_(self%id_l,self%id_b,self%vb*l/(self%k4+min67+l)*(b+self%b0))
   
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine fasham_do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Fasham model
!
! !INTERFACE:
   subroutine fasham_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
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
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_fasham),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                   :: p,z,b,d,n,a,l,par
   REALTYPE                   :: d_p,d_z,d_b,d_d,d_n,d_a,d_l,d_par
   REALTYPE                   :: ff,fac,min67
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_p,p) ! phytoplankton
   _GET_STATE_(self%id_z,z) ! zooplankton
   _GET_STATE_(self%id_b,b) ! bacteria
   _GET_STATE_(self%id_d,d) ! detritus
   _GET_STATE_(self%id_n,n) ! nitrate
   _GET_STATE_(self%id_a,a) ! ammonia
   _GET_STATE_(self%id_l,l) ! labile dissolved organic nitrogen
   
   ! Retrieve local photosynthetically active radiation.
   _GET_DEPENDENCY_(self%id_par,par)

   ! Calculate intermediate quantities that will be used multiple times.
   ff = self%vp*self%alpha*par/sqrt(self%vp**2+self%alpha**2*par**2) 
   if (p.eq._ZERO_ .and. b.eq._ZERO_ .and. d.eq._ZERO_) then
      fac = _ZERO_
   else
      fac = (z+self%z0)/(self%k3*(self%r1*p+self%r2*b+self%r3*d)+  &
                         self%r1*p**2+self%r2*b**2+self%r3*d**2)
   end if
   min67 = min(a,self%eta*l)

   d_p = ff*(n/self%k1+a/self%k2)/(_ONE_+n/self%k1+a/self%k2)*(p*(_ONE_-self%gamma)+self%p0) &
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
   d_d =   (_ONE_-self%beta)*self%gmax*self%r1*p**2*fac &
         + (_ONE_-self%beta)*self%gmax*self%r2*b**2*fac &
         - self%beta        *self%gmax*self%r3*d**2*fac &
         - self%mu4*d                                   &
         + self%mu1*(p+self%p0)/(self%k5+p+self%p0)*p   &
         + (_ONE_-self%epsi-self%delta)*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z
   d_n = - ff*n/self%k1/(_ONE_+n/self%k1+a/self%k2)*(p+self%p0)
   d_a = - ff*a/self%k2/(_ONE_+n/self%k1+a/self%k2)*(p+self%p0) &
         - self%vb*min67/(self%k4+min67+l)*(b+self%b0)          &
         + self%mu3*b                                           &
         + self%epsi*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z
   d_l =   self%gamma*ff*(n/self%k1+a/self%k2)/(_ONE_+n/self%k1+a/self%k2)*p &
         + self%mu4*d                                                        &
         + self%delta*self%mu2*(z+self%z0)/(self%k6+z+self%z0)*z             &
         - self%vb*l/(self%k4+min67+l)*(b+self%b0)

   _SET_ODE_(self%id_p,d_p)
   _SET_ODE_(self%id_z,d_z)
   _SET_ODE_(self%id_d,d_d)
   _SET_ODE_(self%id_b,d_b)
   _SET_ODE_(self%id_n,d_n)
   _SET_ODE_(self%id_a,d_a)
   _SET_ODE_(self%id_l,d_l)
   
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine fasham_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   pure subroutine fasham_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   type (type_fasham), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: p
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_p,p) ! phytoplankton
   
   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(self%p0+p))

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_
   
   end subroutine fasham_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   pure subroutine fasham_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   type (type_fasham), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: p,z,b,d,n,a,l
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_p,p) ! phytoplankton
   _GET_STATE_(self%id_z,z) ! zooplankton
   _GET_STATE_(self%id_b,b) ! bacteria
   _GET_STATE_(self%id_d,d) ! detritus
   _GET_STATE_(self%id_n,n) ! nitrate
   _GET_STATE_(self%id_a,a) ! ammonia
   _GET_STATE_(self%id_l,l) ! labile dissolved organic nitrogen
   
   ! Total nutrient is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totN,p+z+b+d+n+a+l)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine fasham_get_conserved_quantities
!EOC

!-----------------------------------------------------------------------

   end module fabm_fasham

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

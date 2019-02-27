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

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm\_gotm\_ergom --- GOTM biogeochemical model ERGOM
!
! !INTERFACE:
   module gotm_ergom
!
! !DESCRIPTION:
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
!
! !USES:
   use fabm_types

   implicit none

   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_gotm_ergom
!
! !REVISION HISTORY:
!  Author(s):
!
! !PUBLIC_DERIVED_TYPES:
   type,extends(type_base_model) :: type_gotm_ergom
! Variable identifiers
      type (type_state_variable_id)        :: id_p1,id_p2,id_p3,id_zo,id_de,id_am,id_ni,id_po,id_o2
      type (type_bottom_state_variable_id) :: id_fl
      type (type_dependency_id)            :: id_par,id_temp,id_salt
      type (type_horizontal_dependency_id) :: id_I_0,id_wind,id_taub
      type (type_diagnostic_variable_id)   :: id_dPAR,id_GPP,id_NCP,id_PPR,id_NPR
! Model parameters
      real(rk) :: sfl_po,sfl_am,sfl_ni,p10,p20,p30,zo0,kc,i_min,r1max,r2max,r3max,alpha1,alpha2,alpha3,lpa,lpd
      real(rk) :: tf,tbg,beta_bg,g1max,g2max,g3max,lza,lzd,iv,topt,lan,oan,beta_an,lda,tda,beta_da
      real(rk) :: pvel,sr,s1,s2,s3,s4,a0,a1,a2,lds,lsd,tau_crit,lsa,bsa,ph1,ph2
      logical  :: fluff

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: get_light_extinction
      procedure :: do_surface

  end type
!EOP
!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the ergom model
!
! !INTERFACE:
subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!   Here, the ergom namelist is read and the variables exported by the model are registered with FABM
!
! !USES
!   List any modules used by this routine.
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   class (type_gotm_ergom), intent(inout), target :: self
   integer,                intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk)           :: p1_initial=4.5_rk
   real(rk)           :: p2_initial=4.5_rk
   real(rk)           :: p3_initial=4.5_rk
   real(rk)           :: zo_initial=4.5_rk
   real(rk)           :: de_initial=4.5_rk
   real(rk)           :: am_initial=4.5_rk
   real(rk)           :: ni_initial=4.5_rk
   real(rk)           :: po_initial=4.5_rk
   real(rk)           :: o2_initial=4.5_rk
   real(rk)           :: sfl_po=0.0015_rk
   real(rk)           :: sfl_am=0.07_rk
   real(rk)           :: sfl_ni=0.09_rk
   logical            :: fluff=.false.
   real(rk)           :: fl_initial=0.0_rk
   real(rk)           :: p10=0.0225_rk
   real(rk)           :: p20=0.0225_rk
   real(rk)           :: p30=0.0225_rk
   real(rk)           :: zo0=0.0225_rk
   real(rk)           :: w_p1=-1.157407e-05_rk
   real(rk)           :: w_p2=-5.787037e-05_rk
   real(rk)           :: w_p3=-5.787037e-05_rk
   real(rk)           :: w_de=-3.0_rk
   real(rk)           :: kc=0.03_rk
   real(rk)           :: i_min=25.0_rk
   real(rk)           :: r1max=2.0_rk
   real(rk)           :: r2max=0.7_rk
   real(rk)           :: r3max=0.5_rk
   real(rk)           :: alpha1=1.35_rk
   real(rk)           :: alpha2=0.675_rk
   real(rk)           :: alpha3=0.5_rk
   real(rk)           :: lpa=0.01_rk
   real(rk)           :: lpd=0.02_rk
   real(rk)           :: tf=10.0_rk
   real(rk)           :: tbg=14.0_rk
   real(rk)           :: beta_bg=1._rk
   real(rk)           :: g1max=0.5_rk
   real(rk)           :: g2max=0.5_rk
   real(rk)           :: g3max=0.25_rk
   real(rk)           :: lza=0.0666666666_rk
   real(rk)           :: lzd=0.1333333333_rk
   real(rk)           :: iv=0.24444444_rk
   real(rk)           :: topt=20.0_rk
   real(rk)           :: lan=0.1_rk
   real(rk)           :: oan=0.01_rk
   real(rk)           :: beta_an=0.11_rk
   real(rk)           :: lda=0.003_rk
   real(rk)           :: tda=13.0_rk
   real(rk)           :: beta_da=20.0_rk
   real(rk)           :: pvel=5.0_rk
   real(rk)           :: sr=0.0625_rk
   real(rk)           :: s1=5.3_rk
   real(rk)           :: s2=6.625_rk
   real(rk)           :: s3=8.625_rk
   real(rk)           :: s4=2.0_rk
   real(rk)           :: a0=31.25_rk
   real(rk)           :: a1=14.603_rk
   real(rk)           :: a2=0.4025_rk
   real(rk)           :: lds=3.5_rk
   real(rk)           :: lsd=25.0_rk
   real(rk)           :: tau_crit=0.07_rk
   real(rk)           :: lsa=0.001_rk
   real(rk)           :: bsa=0.15_rk
   real(rk)           :: ph1=0.15_rk
   real(rk)           :: ph2=0.1_rk
   real(rk),parameter :: secs_pr_day=86400.0_rk

   namelist /gotm_ergom/ p1_initial,p2_initial,p3_initial,zo_initial,  &
                        de_initial,am_initial,ni_initial,po_initial,  &
                        o2_initial,sfl_po,sfl_am,sfl_ni,fluff,        &
                        fl_initial,p10,p20,p30,zo0,w_p1,w_p2,w_p3,    &
                        w_de,kc,i_min,r1max,r2max,r3max,alpha1,       &
                        alpha2,alpha3,lpa,lpd,tf,tbg,beta_bg,g1max,   &
                        g2max,g3max,lza,lzd,iv,topt,lan,oan,beta_an,  &
                        lda,tda,beta_da,pvel,sr,s1,s2,s3,s4,a0,a1,a2, &
                        lds,lsd,tau_crit,lsa,bsa,ph1,ph2
!EOP
!-----------------------------------------------------------------------
!BOC
!  Read the namelist
   if (configunit>0) read(configunit,nml=gotm_ergom,err=99,end=100)

!  Store parameter values in our own derived type
!  NB! All rates must be provided in values per day,
!  and are converted here to values per second
   self%sfl_po   = sfl_po
   self%sfl_ni   = sfl_ni
   self%sfl_am   = sfl_am
   self%fluff    = fluff
   self%p10      = p10
   self%p20      = p20
   self%p30      = p30
   self%zo0      = zo0
   self%kc       = kc
   self%i_min    = i_min
   self%r1max    = r1max/secs_pr_day
   self%r2max    = r2max/secs_pr_day
   self%r3max    = r3max/secs_pr_day
   self%alpha1   = alpha1
   self%alpha2   = alpha2
   self%alpha3   = alpha3
   self%lpa      = lpa/secs_pr_day
   self%lpd      = lpd/secs_pr_day
   self%tf       = tf
   self%tbg      = tbg
   self%beta_bg  = beta_bg
   self%g1max    = g1max/secs_pr_day
   self%g2max    = g2max/secs_pr_day
   self%g3max    = g3max/secs_pr_day
   self%lza      = lza/secs_pr_day
   self%lzd      = lzd/secs_pr_day
   self%iv       = iv
   self%topt     = topt
   self%lan      = lan/secs_pr_day
   self%oan      = oan
   self%beta_an  = beta_an
   self%lda      = lda/secs_pr_day
   self%tda      = tda
   self%beta_da  = beta_da
   self%pvel     = pvel/secs_pr_day
   self%sr       = sr
   self%s1       = s1
   self%s2       = s2
   self%s3       = s3
   self%s4       = s4
   self%a0       = a0
   self%a1       = a1
   self%a2       = a2
   if (fluff) then
   self%lds      = lds/secs_pr_day
   self%lsd      = lsd/secs_pr_day
   self%tau_crit = tau_crit
   self%lsa      = lsa/secs_pr_day
   self%bsa      = bsa
   self%ph1      = ph1
   self%ph2      = ph2
   end if

!  Register state variables
   call self%register_state_variable(self%id_p1,'dia','mmol n/m**3','diatoms',      &
         p1_initial,minimum=0.0_rk,vertical_movement=w_p1/secs_pr_day)
   call self%register_state_variable(self%id_p2,'fla','mmol n/m**3','flagellates',  &
         p2_initial,minimum=0.0_rk,vertical_movement=w_p2/secs_pr_day)
   call self%register_state_variable(self%id_p3,'cya','mmol n/m**3','cyanobacteria',&
         p3_initial,minimum=0.0_rk,vertical_movement=w_p3/secs_pr_day)
   call self%register_state_variable(self%id_zo,'zoo','mmol n/m**3','zooplankton',  &
         zo_initial,minimum=0.0_rk)
   call self%register_state_variable(self%id_de,'det','mmol n/m**3','detritus',     &
         de_initial,minimum=0.0_rk,vertical_movement=w_de/secs_pr_day)
   call self%register_state_variable(self%id_am,'amm','mmol n/m**3','ammonium',     &
         am_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_ni,'nit','mmol n/m**3','nitrate',      &
         ni_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_po,'pho','mmol p/m**3','phosphate',    &
         po_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_o2,'oxy','mmol o2/m**3','oxygen',o2_initial)
   if (self%fluff) call self%register_state_variable(self%id_fl,'flf','mmol n/m**2','fluff', &
         fl_initial,minimum=0.0_rk)

! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2','photosynthetically active radiation',   &
         time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_GPP,'GPP','mmol/m**3','gross primary production',            &
         time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_NCP,'NCP','mmol/m**3','net community production',            &
         time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_PPR,'PPR','mmol/m**3/d','gross primary production rate',     &
         time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_NPR,'NPR','mmol/m**3/d','net community production rate',     &
         time_treatment=time_treatment_averaged)

! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_wind,standard_variables%wind_speed)
   if (self%fluff) call self%register_dependency(self%id_taub,standard_variables%bottom_stress)

   return

 99 call self%fatal_error('gotm_ergom_init','Error reading namelist gotm_ergom.')

100 call self%fatal_error('gotm_ergom_init','Namelist gotm_ergom was not found.')

   END subroutine initialize
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of ergom model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   class(type_gotm_ergom), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_
!
!
! !LOCAL VARIABLES:
    real(rk)        :: p1,p2,p3,zo,am,ni,po,de,o2,par,I_0
    real(rk)        :: iopt,ppi,temp,psum,llda,llan,lp,r1,r2,r3
    real(rk)        :: wo=30.0_rk,wn=0.1_rk
    real(rk)        :: thopnp,thomnp,thomnm,thsum
    real(rk)        :: secs_pr_day=86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
!   Enter spatial_loops (if any)
    _LOOP_BEGIN_

!   Retrieve current (local) state variable values
    _GET_(self%id_p1,p1) !diatoms
    _GET_(self%id_p2,p2) !flagellates
    _GET_(self%id_p3,p3) !cyanobacteria
    _GET_(self%id_zo,zo) !zooplankton
    _GET_(self%id_de,de) !detritus
    _GET_(self%id_am,am) !ammonium
    _GET_(self%id_ni,ni) !nitrate
    _GET_(self%id_po,po) !phosphate
    _GET_(self%id_o2,o2) !oxygen

!   Retrieve current environmental conditions
    _GET_   (self%id_par,par) ! local photosynthetically active radiation
    _GET_HORIZONTAL_(self%id_I_0,I_0) ! surface short wave radiation
    _GET_   (self%id_temp,temp) !local water temperature

!   Light acclimation formulation based on surface light intensity
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

  _SET_ODE_(self%id_p1,r1-fpz(self,self%g1max,temp,self%topt,psum)*p1/psum*(zo+self%zo0)-lp*p1)
  _SET_ODE_(self%id_p2,r2-fpz(self,self%g2max,temp,self%topt,psum)*p2/psum*(zo+self%zo0)-lp*p2)
  _SET_ODE_(self%id_p3,r3-fpz(self,self%g3max,temp,self%topt,psum)*p3/psum*(zo+self%zo0)-lp*p3)
  _SET_ODE_(self%id_zo,(fpz(self,self%g1max,temp,self%topt,psum)*p1+fpz(self,self%g2max,temp,self%topt,psum)*p2+fpz(self,self%g3max,temp,self%topt,psum)*p3)*(zo+self%zo0)/psum-self%lza*zo*(zo+self%zo0)-self%lzd*zo*(zo+self%zo0))
  _SET_ODE_(self%id_de,self%lpd*p1+self%lpd*p2+self%lpd*p3+self%lzd*zo*(zo+self%zo0)-llda*de)
  _SET_ODE_(self%id_am,llda*de-llan*am-am/(am+ni)*r1-am/(am+ni)*r2+self%lpa*(p1+p2+p3)+self%lza*zo*(zo+self%zo0))
  _SET_ODE_(self%id_ni,llan*am-ni/(am+ni)*r1-ni/(am+ni)*r2-self%s1*llda*de*thomnp)
  _SET_ODE_(self%id_po,self%sr*(-r1-r2-r3+llda*de+self%lpa*(p1+p2+p3)+self%lza*zo*(zo+self%zo0)))
  _SET_ODE_(self%id_o2,self%s2*(am/(am+ni)*(r1+r2)+r3)-self%s4*(llan*am)+self%s3*(ni/(am+ni)*(r1+r2))-self%s2*(thopnp+thomnm)*llda*de-self%s2*(self%lpa*(p1+p2+p3)+self%lza*zo*(zo+self%zo0)))

! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP ,r1+r2+r3)
   _SET_DIAGNOSTIC_(self%id_NCP ,r1+r2+r3 - self%lpa*(p1+p2+p3))
   _SET_DIAGNOSTIC_(self%id_PPR ,(r1+r2+r3)*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR ,(r1+r2+r3 - self%lpa*(p1+p2+p3))*secs_pr_day)

!   Leave spatial loops (if any)
   _LOOP_END_

   END subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !DESCRIPTION:

! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_gotm_ergom), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                     :: p1,p2,p3,de
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p1,p1) ! diatoms
   _GET_(self%id_p2,p2) !flagellates
   _GET_(self%id_p3,p3) !cyanobacteria
   _GET_(self%id_de,de) ! detritus


   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(self%p10+self%p20+self%p30+p1+p2+p3+de))

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
! Calculating the benthic fluxes
!
   implicit none

! !INPUT PARAMETERS:
   class (type_gotm_ergom), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk)                   :: fl,amb,nib,pob,deb,oxb,taub,temp
   real(rk)                   :: llds,llsd,llsa,wo=30.0_rk,wn=0.1_rk,dot2=0.2_rk
   real(rk)                   :: thopnp,thomnp,thomnm,thsum

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   if (self%fluff) then
      _GET_(self%id_am,amb)
      _GET_(self%id_de,deb)
      _GET_(self%id_ni,nib)
      _GET_(self%id_po,pob)
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

        if (self%tau_crit .gt. taub) then
            llds=self%lds*(self%tau_crit-taub)/self%tau_crit
         else
            llds=0.
         end if
         if (self%tau_crit .lt. taub) then
            llsd=self%lsd*(taub-self%tau_crit)/self%tau_crit
         else
            llsd=0.0_rk
         end if

      _SET_BOTTOM_ODE_(self%id_fl,llds*deb-llsd*fl-llsa*fl-th(oxb,wo,0.0_rk,1.0_rk)*llsa*fl)
      _SET_BOTTOM_EXCHANGE_(self%id_de,-llds*deb+llsd*fl)
      _SET_BOTTOM_EXCHANGE_(self%id_am,llsa*fl)
      _SET_BOTTOM_EXCHANGE_(self%id_ni,-self%s1*thomnp*llsa*fl)
      _SET_BOTTOM_EXCHANGE_(self%id_po,self%sr*(1.0_rk-self%ph1*th(oxb,wo,0.0_rk,1.0_rk)*yy(self%ph2,oxb))*llsa*fl)
      _SET_BOTTOM_EXCHANGE_(self%id_o2,-(self%s4+self%s2*(thopnp+thomnm))*llsa*fl)
   end if

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
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
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!  type (type_gotm_ergom), intent(in) :: self
  real(rk), intent(in)                 :: t,s
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  real(rk)                 :: tk
  real(rk)                 :: aa1=-173.4292_rk
  real(rk)                 :: aa2=249.6339_rk
  real(rk)                 :: a3=143.3483_rk
  real(rk)                 :: a4=-21.8492_rk
  real(rk)                 :: b1=-0.033096_rk
  real(rk)                 :: b2=0.014259_rk
  real(rk)                 :: b3=-0.001700_rk
  real(rk)                 :: kelvin=273.16_rk
  real(rk)                 :: mol_per_liter=44.661_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   tk=(t+kelvin)*0.01_rk
   osat_weiss=exp(aa1+aa2/tk+a3*log(tk)+a4*tk    &
              +s*(b1+(b2+b3*tk)*tk))*mol_per_liter
   return
   end function osat_weiss
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the ergom model
!
! !INTERFACE:

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_gotm_ergom),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

  real(rk)             :: temp,wnd,salt,o2,ni,am,po
  real(rk),parameter           :: secs_pr_day=86400.0_rk
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
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
   _GET_(self%id_po,po)

!  Calculation of the surface oxygen flux
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
   _SET_SURFACE_EXCHANGE_(self%id_o2,flo2)
   else
      flo2 = self%pvel*(self%a0*(self%a1-self%a2*temp)-o2)
   _SET_SURFACE_EXCHANGE_(self%id_o2,flo2)
   end if

   _SET_SURFACE_EXCHANGE_(self%id_ni,self%sfl_ni/secs_pr_day)
   _SET_SURFACE_EXCHANGE_(self%id_am,self%sfl_am/secs_pr_day)
   _SET_SURFACE_EXCHANGE_(self%id_po,self%sfl_po/secs_pr_day)

   _HORIZONTAL_LOOP_END_
   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Step function
!
! !INTERFACE:
   pure real(rk) function th(x,w,min,max)
!
! !DESCRIPTION:
! Instead of the
! heavyside switches used by \cite{Neumannetal2002}, we apply here a smoothed
! {\it tangens hyperbolicus} transition with prescribed width $x_w$:
! \begin{equation}\label{theta}
! \theta (x,x_w,y_{\min},y_{\max})= y_{\min}+(y_{\max}-y_{\min})
! \frac12\left(1-\tanh \left(\frac{x}{x_w}   \right)      \right).
! \end{equation}
!
! !USES
!   List any modules used by this routine.
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
   ! type(type_gotm_ergom), INTENT(IN) :: self
    real(rk), intent(in)            :: x,w,min,max
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
    if (w .gt. 1.e-10_rk) then
       th=min+(max-min)*0.5_rk*(1.0_rk+tanh(x/w))
    else
       if(x .gt. 0.0_rk) then
          th=1.0_rk
       else
          th=0.0_rk
       end if
    end if
   RETURN
   END function th
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:
   pure real(rk) function yy(a,x)
!
! !DESCRIPTION:
! This is a squared Michaelis-Menten type of limiter:
! \begin{equation}\label{Y}
! Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
! \end{equation}
!
! !USES:
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  real(rk),intent(in)        :: a,x
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC
  yy=x**2/(a**2+x**2)
  RETURN
  END function yy
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
  pure real(rk) function fpz(self,g,t,topt,psum)
!
! !DESCRIPTION:
! The Ivlev formulation for zooplankton grazing on the three phytoplankton
! classes $c_1$, $c_2$, and $c_3$ is given here as a function:
! \begin{equation}\label{neu_di4}
! d_{i,4}=g_i^{\max}\left(1+\frac{T^2}{T_{opt}^2}\exp
! \left(1-\frac{2T}{T_{opt}} \right)\right)
! \left( 1-\exp\left(-I_v^2 \left( \sum_{
! j=1}^3c_j \right)^2\right)  \right)
! \frac{c_i}{\sum_{j=1}^3c_j}\left( c_4+c_4^{\min} \right)
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   class (type_gotm_ergom), intent(in) :: self
   real(rk), intent(in)                :: g,t,topt,psum
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fpz=g*(1.0_rk+t**2/topt**2*exp(1.0_rk-2.0_rk*t/topt))*               &
        (1.0_rk-exp(-self%iv**2*psum**2))
   return
   end function fpz
!EOC

!-----------------------------------------------------------------------

  END MODULE gotm_ergom

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

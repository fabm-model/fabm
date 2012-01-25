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
! !MODULE: fabm\_ergom --- IOW biogeochemical model ERGOM 
!
! !INTERFACE:
   module fabm_iow_ergom
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
! model is given in section \ref{sec:bio-iow-details}.
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
   use fabm_driver

   implicit none

   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_iow_ergom, iow_ergom_init, iow_ergom_do, iow_ergom_get_light_extinction, &
   iow_ergom_get_surface_exchange,iow_ergom_do_benthos
!
! !PUBLIC DATA MEMBERS:
!
! !REVISION HISTORY:
!  Author(s):
!
! !PUBLIC_DERIVED_TYPES:
  type type_iow_ergom
! Variable identifiers
      _TYPE_STATE_VARIABLE_ID_      :: id_p1,id_p2,id_p3,id_zo,id_de,id_am,id_ni,id_po,id_o2,id_fl
      _TYPE_DEPENDENCY_ID_          :: id_par,id_I_0,id_temp,id_salt,id_wind,id_taub
      _TYPE_DIAGNOSTIC_VARIABLE_ID_ :: id_dPAR,id_GPP,id_NCP,id_PPR,id_NPR
! Model parameters
      REALTYPE :: sfl_po,sfl_am,sfl_ni,p10,p20,p30,zo0,kc,i_min,r1max,r2max,r3max,alpha1,alpha2,alpha3,lpa,lpd
      REALTYPE :: tf,tbg,beta_bg,g1max,g2max,g3max,lza,lzd,iv,topt,lan,oan,beta_an,lda,tda,beta_da
      REALTYPE :: pvel,sr,s1,s2,s3,s4,a0,a1,a2,lds,lsd,tau_crit,lsa,bsa,ph1,ph2
      logical  :: fluff
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
   subroutine iow_ergom_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!   Here, the ergom namelist is read and the variables exported by the model are registered with FABM
!
! !USES
!   List any modules used by this routine.
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type(type_iow_ergom),intent(out)     ::self
   type(type_model_info),intent(inout) :: modelinfo
   integer,intent(in)                  ::namlst
!
! !LOCAL VARIABLES:
   REALTYPE           :: p1_initial=4.5
   REALTYPE           :: p2_initial=4.5
   REALTYPE           :: p3_initial=4.5
   REALTYPE           :: zo_initial=4.5
   REALTYPE           :: de_initial=4.5
   REALTYPE           :: am_initial=4.5
   REALTYPE           :: ni_initial=4.5
   REALTYPE           :: po_initial=4.5
   REALTYPE           :: o2_initial=4.5
   REALTYPE           :: sfl_po=0.0015
   REALTYPE           :: sfl_am=0.07
   REALTYPE           :: sfl_ni=0.09
   logical            :: fluff=.false.
   REALTYPE           :: fl_initial=0.0
   REALTYPE           :: p10=0.0225
   REALTYPE           :: p20=0.0225
   REALTYPE           :: p30=0.0225
   REALTYPE           :: zo0=0.0225
   REALTYPE           :: w_p1=-1.157407e-05
   REALTYPE           :: w_p2=-5.787037e-05
   REALTYPE           :: w_p3=-5.787037e-05
   REALTYPE           :: w_de=-3.
   REALTYPE           :: kc=0.03
   REALTYPE           :: i_min=25.
   REALTYPE           :: r1max=2.0
   REALTYPE           :: r2max=0.7
   REALTYPE           :: r3max=0.5
   REALTYPE           :: alpha1=1.35
   REALTYPE           :: alpha2=0.675
   REALTYPE           :: alpha3=0.5
   REALTYPE           :: lpa=0.01
   REALTYPE           :: lpd=0.02
   REALTYPE           :: tf=10.
   REALTYPE           :: tbg=14.
   REALTYPE           :: beta_bg=1.
   REALTYPE           :: g1max=0.5
   REALTYPE           :: g2max=0.5
   REALTYPE           :: g3max=0.25
   REALTYPE           :: lza=0.0666666666
   REALTYPE           :: lzd=0.1333333333
   REALTYPE           :: iv=0.24444444
   REALTYPE           :: topt=20.
   REALTYPE           :: lan=0.1
   REALTYPE           :: oan=0.01
   REALTYPE           :: beta_an=0.11
   REALTYPE           :: lda=0.003
   REALTYPE           :: tda=13.
   REALTYPE           :: beta_da=20.
   REALTYPE           :: pvel=5.
   REALTYPE           :: sr=0.0625
   REALTYPE           :: s1=5.3
   REALTYPE           :: s2=6.625
   REALTYPE           :: s3=8.625
   REALTYPE           :: s4=2.0
   REALTYPE           :: a0=31.25
   REALTYPE           :: a1=14.603
   REALTYPE           :: a2=0.4025
   REALTYPE           :: lds=3.5
   REALTYPE           :: lsd=25.0
   REALTYPE           :: tau_crit=0.07
   REALTYPE           :: lsa=0.001
   REALTYPE           :: bsa=0.15
   REALTYPE           :: ph1=0.15
   REALTYPE           :: ph2=0.1

   REALTYPE,parameter           :: secs_pr_day=86400.
   namelist /iow_ergom/ p1_initial,p2_initial,p3_initial,zo_initial,  &
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
   read(namlst,nml=iow_ergom,err=99,end=100)

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
   self%id_p1=register_state_variable(modelinfo,'dia','mmol n/m**3','diatoms',p1_initial,minimum=_ZERO_,vertical_movement=w_p1/secs_pr_day)
   self%id_p2=register_state_variable(modelinfo,'fla','mmol n/m**3','flagellates',p2_initial,minimum=_ZERO_,vertical_movement=w_p2/secs_pr_day)
   self%id_p3=register_state_variable(modelinfo,'cya','mmol n/m**3','cyanobacteria',p3_initial,minimum=_ZERO_,vertical_movement=w_p3/secs_pr_day)
   self%id_zo=register_state_variable(modelinfo,'zoo','mmol n/m**3','zooplankton',zo_initial,minimum=_ZERO_)
   self%id_de=register_state_variable(modelinfo,'det','mmol n/m**3','detritus',de_initial,minimum=_ZERO_,vertical_movement=w_de/secs_pr_day)
   self%id_am=register_state_variable(modelinfo,'amm','mmol n/m**3','ammonium',am_initial,minimum=_ZERO_,no_river_dilution=.true.)
   self%id_ni=register_state_variable(modelinfo,'nit','mmol n/m**3','nitrate',ni_initial,minimum=_ZERO_,no_river_dilution=.true.)
   self%id_po=register_state_variable(modelinfo,'pho','mmol p/m**3','phosphate',po_initial,minimum=_ZERO_,no_river_dilution=.true.)
   self%id_o2=register_state_variable(modelinfo,'oxy','mmol o2/m**3','oxygen',o2_initial)
   if (self%fluff) self%id_fl=register_state_variable(modelinfo,'flf','mmol n/m**2','fluff',fl_initial,benthic=.true.,minimum=_ZERO_)

! Register diagnostic variables
   self%id_dPAR = register_diagnostic_variable(modelinfo,'PAR','W/m**2','photosynthetically active radiation',time_treatment=time_treatment_averaged)
   self%id_GPP  = register_diagnostic_variable(modelinfo,'GPP','mmol/m**3','gross primary production',time_treatment=time_treatment_step_integrated)
   self%id_NCP  = register_diagnostic_variable(modelinfo,'NCP','mmol/m**3','net community production',time_treatment=time_treatment_step_integrated)
   self%id_PPR  = register_diagnostic_variable(modelinfo,'PPR','mmol/m**3/d','gross primary production rate',time_treatment=time_treatment_averaged)
   self%id_NPR  = register_diagnostic_variable(modelinfo,'NPR','mmol/m**3/d','net community production rate',time_treatment=time_treatment_averaged)

! Register environmental dependencies
   self%id_par = register_dependency(modelinfo, varname_par)
   self%id_I_0 = register_dependency(modelinfo, varname_par_sf, shape=shape_hz)
   self%id_temp = register_dependency(modelinfo, varname_temp)
   self%id_salt = register_dependency(modelinfo, varname_salt)
   self%id_wind = register_dependency(modelinfo, varname_wind_sf,shape=shape_hz)
   if (self%fluff) self%id_taub=register_dependency(modelinfo,varname_taub,shape=shape_hz)

   RETURN

 99 call fatal_error('iow_ergom_init','Error reading namelist iow_ergom.')

100 call fatal_error('iow_ergom_init','Namelist iow_ergom was not found.')

   END subroutine iow_ergom_init
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of ergom model
!
! !INTERFACE:
   subroutine iow_ergom_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
! The right hand sides of the \cite{Neumannetal2002} biogeochemical model are
! coded in this soubroutine.
! First of all, based on (\ref{theta}) and (\ref{Y}),
! we construct limiters for chemical
! reactions which depend on the availability of oxygen ($c_9$) and nitrate
! ($c_7$) and have to add up to unity:
! \begin{equation}\label{limits}
! \begin{array}{rcl}
! l^+_+ &=& \theta(c_9,c_9^t,0,1)Y(c_7^t,c_7), \\ \\
! l^-_+ &=& \theta(-c_9,c_9^t,0,1)Y(c_7^t,c_7), \\ \\
! l^-_- &=& \theta(-c_9,c_9^t,0,1)(1-Y(c_7^t,c_7)), \\ \\
! L^+_+ &=& \frac{l^+_+}{l^+_+ + l^-_+ + l^-_-}, \\ \\
! L^-_+ &=& \frac{l^-_+}{l^+_+ + l^-_+ + l^-_-}, \\ \\
! L^-_- &=& \frac{l^-_-}{l^+_+ + l^-_+ + l^-_-}. \\ \\
! \end{array}
! \end{equation}
!
! Mortality of the three phytoplankton classes $c_i$, $i=1,\dots,3$:
! \begin{equation}\label{neu_di5}
! d_{i,5}=l_{PD} c_i
! \end{equation}
!
! Respiration of the three phytoplankton classes $c_i$, $i=1,\dots,3$
! into ammonium:
! \begin{equation}\label{neu_di6}
! d_{i,6}=l_{PA} c_i
! \end{equation}
!
! Zooplankton mortality:
! \begin{equation}\label{neu_d45}
! d_{4,5}=l_{ZD}(c_4+c_4^{\min})c_4
! \end{equation}
!
! Zooplankton exudation into ammonium:
! \begin{equation}\label{neu_d46}
! d_{4,6}=l_{ZA}(c_4+c_4^{\min})c_4
! \end{equation}
!
! Detritus mineralisation:
! \begin{equation}\label{neu_d56}
! d_{5,6}=L_{DA}c_5
! \end{equation}
! with
! \begin{equation}\label{LDA}
! L_{DA} = l_{DA} \left(1+\beta_{DA}Y(T_{DA},T)\right).
! \end{equation}
!
! Ammonium uptake by phytoplankta $c_i$, $i=1,2$:
! \begin{equation}\label{neu_d6i}
! d_{6,i}=R_i\frac{c_6}{c_6+c_7}\left(c_i+c_i^{\min} \right)
! \end{equation}
! with the growth rate for diatoms,
! \begin{equation}\label{r1}
! R_1=r_1^{\max} \min\left\{
! Y(\alpha_1,c_6+c_7), Y(s_R\alpha_1,c_8), PPI
! \right\}
! \end{equation}
! and the growth rate for flagellates,
! \begin{equation}\label{r2}
! R_2=r_2^{\max}\left(1+Y\left(T_f,T \right)\right)\min\left\{
! Y(\alpha_2,c_6+c_7),Y(s_R\alpha_2,c_8),PPI
! \right\}.
! \end{equation}
! Here,
! \begin{equation}\label{ppi}
! PPI=\frac{I_{PAR}}{I_{opt}}\exp\left(1-\frac{I_{PAR}}{I_{opt}}  \right)
! \end{equation}
! with
! \begin{equation}\label{iopt}
! I_{opt}=\max\left\{\frac{I_0}{4},I_{\min}   \right\}
! \end{equation}
! and $I_{PAR}$ from (\ref{light}).
!
! Nitrification of ammonium to nitrate:
! \begin{equation}\label{neu_d67}
! d_{6,7}=L_{AN}c_6
! \end{equation}
! with
! \begin{equation}\label{LAN}
! L_{AN}=l_{AN}\theta(c_9,0,0,1)\frac{c_9}{O_{AN}+c_9}\exp\left(\beta_{AN}T\right).
! \end{equation}
!
! Nitrate uptake by phytoplankta $c_i$, $i=1,2$:
! \begin{equation}\label{neu_d7i}
! d_{7,i}=R_i\frac{c_7}{c_6+c_7}\left(c_i+c_i^{\min} \right).
! \end{equation}
!
! Settling of detritus into sediment:
! \begin{equation}\label{neu_d510}
! d_{5,10}=l_{DS} \frac{c_5}{h_1}\delta_{k,1}
! \end{equation}
!
! Mineralisation of sediment into ammonium:
! \begin{equation}\label{neu_d106}
! d_{10,6}=L_{SA} c_{10}
! \end{equation}
! with
! \begin{equation}\label{LSA}
! L_{SA}=l_{SA} \exp\left(\beta_{SA}T \right) \theta(c_9,c_9^t,0.2,1)
! \end{equation}
!
! From the above sink terms, respective source terms are calculated by means of (\ref{eq:am:symmetry}),
! except for settling of detritus into sediment and mineralisation of sediment into
! ammonium, for which we have:
! \begin{equation}\label{neu_p105}
! p_{10,5}=h_1 d_{5,10}, \quad p_{6,10}=\frac{d_{10,6}}{h_1}.
! \end{equation}
!
! Denitrification in water column:
! \begin{equation}\label{neu_d77}
! d_{7,7}=s_1 \left(L_{DA} c_5 +L_{SA} \frac{c_{10}}{h_1}\delta_{k,1} \right)L^-_+.
! \end{equation}
!
! Denitrification in sediment:
! \begin{equation}\label{neu_d1010}
! d_{10,10}=\theta(c_9,c_9^t,0,1) L_{SA} c_{10}
! \end{equation}
!
! Phosphorus uptake by the three phytoplankton classes $c_i$, $i=1,\dots,3$:
! \begin{equation}\label{neu_d88}
! d_{8,8}=s_R \left(\sum_{j=1}^3 R_j \left(c_j+c_j^{\min}   \right)     \right).
! \end{equation}
!
! Nitrogen fixation:
! \begin{equation}\label{neu_p33}
! p_{3,3}=R_3\left(c_3+c_3^{\min}\right)
! \end{equation}
! with
! \begin{equation}\label{r3}
! R_3=r_3^{\max}\frac{1}{1+\exp\left(\beta_{bg}\left(T_{bg}-T  \right)   \right)}\min\left\{Y\left(s_R\alpha_3,c_8\right),PPI   \right\}
! \end{equation}
!
! Respiration of the three phytoplankton classes $c_i$, $i=1,\dots,3$
! into phosphorus:
! \begin{equation}\label{neu_p8i}
! p_{8,i}=s_R d_{i,6}.
! \end{equation}
!
! Zooplankton exudation into phosphorus:
! \begin{equation}\label{neu_p84}
! p_{8,4}=s_R d_{4,6}.
! \end{equation}
!
! Phosphate release due to detritus mineralisation:
! \begin{equation}\label{neu_p85}
! p_{8,5}=s_R d_{5,6}.
! \end{equation}
!
! Oxygen production due to ammonium uptake by phytoplankton classes $c_i$, $i=1,2$and nitrification of ammonium into nitrate:
! \begin{equation}\label{neu_p96}
! p_{9,6}= s_2 \left(d_{6,1}+d_{6,2} \right)-s_4 d_{6,7}.
! \end{equation}
!
! Oxygen production due to nitrate uptake by phytoplankton classes $c_i$, $i=1,2$:
! \begin{equation}\label{neu_p97}
! p_{9,7}= s_3 \left(d_{7,1}+d_{7,2} \right).
! \end{equation}
!
! Oxygen production due to nitrogen fixation by blue-greens:
! \begin{equation}\label{neu_p99}
! p_{9,9}=s_2 p{3,3}
! \end{equation}
!
! Oxygen demand due to
! respiration of the three phytoplankton classes $c_i$, $i=1,\dots,3$:
! \begin{equation}\label{neu_p9i}
! p_{9,i}=-s_2 d_{i,6}.
! \end{equation}
!
! Oxygen demand of zooplankton exudation:
! \begin{equation}\label{neu_p94}
! p_{9,4}=-s_2 d_{4,6}.
! \end{equation}
!
! Oxygen demand of mineralisation of detritus into ammonium:
! \begin{equation}\label{neu_p95}
! p_{9,5}=-s_2\left(L_+^++L_-^- \right) d_{5,6}.
! \end{equation}
!
! Oxygen demand of mineralisation of sediment into ammonium:
! \begin{equation}\label{neu_p910}
! p_{9,10}=-\left(s_4+s_2\left(L_+^++L_-^- \right)\right) \frac{d_{10,6}}{h_1}\delta_{k,1}.
! \end{equation}
!
! Phosphate release due to mineralisation of sediment into ammonium:
! \begin{equation}\label{neu_p88}
! p_{8,8}=s_R(1-p_1\theta\left(c_9,c_9^t,0,1\right)Y(p_2,c_9))\frac{d_{10,6}}{h_1}\delta_{k,1}.
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type(type_iow_ergom), INTENT(IN) :: self
  _DECLARE_FABM_ARGS_DO_RHS_
!
!
! !LOCAL VARIABLES:
    REALTYPE        :: p1,p2,p3,zo,am,ni,po,de,o2,par,I_0
    REALTYPE        :: iopt,ppi,temp,psum,llda,llan,lp,r1,r2,r3
    REALTYPE        :: wo=30.,wn=0.1
    REALTYPE        :: thopnp,thomnp,thomnm,thsum
    REALTYPE        :: secs_pr_day=86400.
!EOP
!-----------------------------------------------------------------------
!BOC
!   Enter spatial_loops (if any)
    _FABM_LOOP_BEGIN_

!   Retrieve current (local) state variable values
    _GET_STATE_(self%id_p1,p1) !diatoms
    _GET_STATE_(self%id_p2,p2) !flagellates
    _GET_STATE_(self%id_p3,p3) !cyanobacteria
    _GET_STATE_(self%id_zo,zo) !zooplankton
    _GET_STATE_(self%id_de,de) !detritus
    _GET_STATE_(self%id_am,am) !ammonium
    _GET_STATE_(self%id_ni,ni) !nitrate
    _GET_STATE_(self%id_po,po) !phosphate
    _GET_STATE_(self%id_o2,o2) !oxygen

!   Retrieve current environmental conditions
    _GET_DEPENDENCY_   (self%id_par,par) ! local photosynthetically active radiation
    _GET_DEPENDENCY_HZ_(self%id_I_0,I_0) ! surface short wave radiation
    _GET_DEPENDENCY_   (self%id_temp,temp) !local water temperature

!   Light acclimation formulation based on surface light intensity
    iopt = max(0.25*I_0,self%i_min)
    ppi = par/iopt*exp(_ONE_-par/iopt)

    thopnp=th( o2,wo,_ZERO_,_ONE_)*yy(wn,ni)
    thomnp=th(-o2,wo,_ZERO_,_ONE_)*yy(wn,ni)
    thomnm=th(-o2,wo,_ZERO_,_ONE_)*(_ONE_-yy(wn,ni))
    thsum=thopnp+thomnp+thomnm
    thopnp=thopnp/thsum
    thomnp=thomnp/thsum
    thomnm=thomnm/thsum

    psum=p1+p2+p3+self%p10+self%p20+self%p30

    llda=self%lda*(_ONE_+self%beta_da*yy(self%tda,temp))
    llan=th(o2,_ZERO_,_ZERO_,_ONE_)*o2/(self%oan+o2)*self%lan*exp(self%beta_an*temp)

    lp=self%lpa+self%lpd

    r1=self%r1max*min(yy(self%alpha1,am+ni),yy(self%sr*self%alpha1,po),ppi)*(p1+self%p10)
    r2=self%r2max*(_ONE_+yy(self%tf,temp))*min(yy(self%alpha2,am+ni),   &
       yy(self%sr*self%alpha2,po),ppi)*(p2+self%p20)
    r3=self%r3max*_ONE_/(_ONE_+exp(self%beta_bg*(self%tbg-temp)))       &
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
   _SET_DIAG_(self%id_dPAR,par)
   _SET_DIAG_(self%id_GPP ,r1/(p1+self%p10)+r2/(p2+self%p20)+r3/(p2+self%p30))
   _SET_DIAG_(self%id_NCP ,r1/(p1+self%p10)+r2/(p2+self%p20)+r3/(p2+self%p30) - self%lpa*(p1+p2+p3))
   _SET_DIAG_(self%id_PPR ,(r1/(p1+self%p10)+r2/(p2+self%p20)+r3/(p2+self%p30))*secs_pr_day)
   _SET_DIAG_(self%id_NPR ,(r1/(p1+self%p10)+r2/(p2+self%p20)+r3/(p2+self%p30) - self%lpa*(p1+p2+p3))*secs_pr_day)

!   Leave spatial loops (if any)
   _FABM_LOOP_END_

   END subroutine iow_ergom_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !DESCRIPTION:

! !INTERFACE:
   pure subroutine iow_ergom_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   type (type_iow_ergom), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: p1,p2,p3,de
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_p1,p1) ! diatoms
   _GET_STATE_(self%id_p2,p2) !flagellates
   _GET_STATE_(self%id_p3,p3) !cyanobacteria
   _GET_STATE_(self%id_de,de) ! detritus


   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(self%p10+self%p20+self%p30+p1+p2+p3+de))

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine iow_ergom_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of benthic_predator model
!
! !INTERFACE:
   subroutine iow_ergom_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!
! !DESCRIPTION:
! This routine calculates the benthic sink and source terms, as well as
! (matching) bottom fluxes for pelagic variables. Both have units mmol/m**2/s.
! Benthic processes are explained in the description of iow\_ergom\_do 
! subroutine.
!
! !INPUT PARAMETERS:
   type (type_iow_ergom),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   REALTYPE                   :: fl,amb,nib,pob,deb,oxb,taub,temp
   !logical                   :: fluff,ltaub
   REALTYPE                   :: llds,llsd,llsa,wo=30.,wn=0.1,dot2=0.2
   REALTYPE                   :: thopnp,thomnp,thomnm,thsum

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _FABM_HZ_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   if (self%fluff) then
   _GET_STATE_(self%id_am,amb)
   _GET_STATE_(self%id_de,deb)
   _GET_STATE_(self%id_ni,nib)
   _GET_STATE_(self%id_po,pob)
   _GET_STATE_(self%id_o2,oxb)
   _GET_STATE_BEN_(self%id_fl,fl)

   _GET_DEPENDENCY_HZ_(self%id_taub,taub)
   _GET_DEPENDENCY_(self%id_temp,temp)

    thopnp=th( oxb,wo,_ZERO_,_ONE_)*yy(wn,nib)
    thomnp=th(-oxb,wo,_ZERO_,_ONE_)*yy(wn,nib)
    thomnm=th(-oxb,wo,_ZERO_,_ONE_)*(_ONE_-yy(wn,nib))
    thsum=thopnp+thomnp+thomnm
    thopnp=thopnp/thsum
    thomnp=thomnp/thsum
    thomnm=thomnm/thsum

  llsa=self%lsa*exp(self%bsa*temp)*(th(oxb,wo,dot2,_ONE_))

  !ltaub=taub**2*1000.

        if (self%tau_crit .gt. taub) then
            llds=self%lds*(self%tau_crit-taub)/self%tau_crit
         else
            llds=0.
         end if
         if (self%tau_crit .lt. taub) then
            llsd=self%lsd*(taub-self%tau_crit)/self%tau_crit
         else
            llsd=_ZERO_
         end if

   _SET_ODE_BEN_(self%id_fl,llds*deb-llsd*fl-llsa*fl-th(oxb,wo,_ZERO_,_ONE_)*llsa*fl)
   _SET_BOTTOM_EXCHANGE_(self%id_de,-llds*deb+llsd*fl)
   _SET_BOTTOM_EXCHANGE_(self%id_am,llsa*fl)
   _SET_BOTTOM_EXCHANGE_(self%id_ni,self%s1*thomnp*llsa*fl)
   _SET_BOTTOM_EXCHANGE_(self%id_po,self%sr*(_ONE_-self%ph1*th(oxb,wo,_ZERO_,_ONE_)*yy(self%ph2,oxb))*llsa*fl)
   _SET_BOTTOM_EXCHANGE_(self%id_o2,-(self%s4+self%s2*(thopnp+thomnm))*llsa*fl)
   end if

   ! Leave spatial loops over the horizontal domain (if any).
   _FABM_HZ_LOOP_END_

   end subroutine iow_ergom_do_benthos
!EOC

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Weiss formula for the saturation oxygen (osat)
!
! !INTERFACE:
   REALTYPE function osat_weiss(t,s)
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
!  type (type_iow_ergom), intent(in) :: self
  REALTYPE, intent(in)                 :: t,s
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE                 :: tk
  REALTYPE                 :: aa1=-173.4292
  REALTYPE                 :: aa2=249.6339
  REALTYPE                 :: a3=143.3483
  REALTYPE                 :: a4=-21.8492
  REALTYPE                 :: b1=-0.033096
  REALTYPE                 :: b2=0.014259
  REALTYPE                 :: b3=-0.001700
  REALTYPE                 :: kelvin=273.16
  REALTYPE                 :: mol_per_liter=44.661
!EOP
!-----------------------------------------------------------------------
!BOC
   tk=(t+kelvin)*0.01
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
  subroutine iow_ergom_get_surface_exchange(self,_FABM_ARGS_GET_SURFACE_EXCHANGE_)
!
! !DESCRIPTION:
! Here, those surface fluxes which have been read from a file are transformed
! to SI units. The surface oxygen flux is calculated by means of the
! following formula:
! \begin{equation}\label{o2flux}
! F^s_9 = p_{vel} \left(O_{sat}-c_9 \right)
! \end{equation}
! The piston velocity, $p_{vel}$, is calculated by using the model of
! \cite{LissMerlivat86}, which includes three regimes (smooth surface,
! rough surface and breaking waves) depending on the magnitude of wind speed,
! $W$:
! \begin{equation}
! \label{p_vel}
! \begin{array}{l}
! p_{vel}=
! \left\{
! \begin{array}{ll}
! 1.003 W / Sc^{0.66} & \mbox{ for } W < 3.6 \mbox{m/s}, \\
! 5.9(2.85 W - 9.65)/Sc^{0.5} & \mbox{ for } 3.6\mbox{m/s} \leq W \leq 13\mbox{m/s}, \\
! 5.9(5.9 W - 49.3) / Sc^{0.5} & \mbox{ for } W > 13 \mbox{m/s}.
! \end{array}
! \right.
! \end{array}
! \end{equation}
! The Schmidt number $Sc$ is defined as ratio between the kinematic viscosity
! and the molecular diffusivity of oxygen. The following expression for $Sc$
! (\cite{Stigebrandt1991}) is applied:
! \begin{equation}
! \label{Sc}
! Sc{=}1450-71T+1.1T^2.
! \end{equation}
! The Weiss formula for the saturation oxygen is used for calculating
! the air-sea surface flux:
! \begin{equation}\label{osat}
! O_{sat}= a_0\left(a_1-a_2T  \right).
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  type(type_iow_ergom),intent(in)       ::self

  _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_

  REALTYPE             :: temp,wnd,salt,o2,ni,am,po
  REALTYPE,parameter           :: secs_pr_day=86400.
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
  REALTYPE                 :: p_vel,sc,flo2
  integer,parameter        :: newflux=1
!EOP
!-----------------------------------------------------------------------
!BOC
   _FABM_HZ_LOOP_BEGIN_

   _GET_DEPENDENCY_(self%id_temp,temp)
   _GET_DEPENDENCY_(self%id_salt,salt)
   _GET_DEPENDENCY_HZ_(self%id_wind,wnd)

   _GET_STATE_(self%id_o2,o2)
   _GET_STATE_(self%id_ni,ni)
   _GET_STATE_(self%id_am,am)
   _GET_STATE_(self%id_po,po)

!  Calculation of the surface oxygen flux
   if (newflux .eq. 1) then
      sc=1450.+(1.1*temp-71.)*temp
      if (wnd .gt. 13.) then
         p_vel = 5.9*(5.9*wnd-49.3)/sqrt(sc)
      else
         if (wnd .lt. 3.6) then
            p_vel = 1.003*wnd/(sc)**(0.66)
         else
            p_vel = 5.9*(2.85*wnd-9.65)/sqrt(sc)
         end if
      end if
      if (p_vel .lt. 0.05) then
          p_vel = 0.05
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

   _FABM_HZ_LOOP_END_
   end subroutine iow_ergom_get_surface_exchange
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Step function
!
! !INTERFACE:
   pure REALTYPE function th(x,w,min,max)
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
   ! type(type_iow_ergom), INTENT(IN) :: self
    REALTYPE, intent(in)            :: x,w,min,max
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
    if (w .gt. 1.e-10) then
       th=min+(max-min)*0.5*(1.+tanh(x/w))
    else
       if(x .gt. _ZERO_) then
          th=_ONE_
       else
          th=_ZERO_
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
   pure REALTYPE function yy(a,x)
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
  REALTYPE,intent(in)        :: a,x
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
  pure REALTYPE function fpz(self,g,t,topt,psum)
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
   type (type_iow_ergom), intent(in) :: self
   REALTYPE, intent(in)                :: g,t,topt,psum
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fpz=g*(_ONE_+t**2/topt**2*exp(_ONE_-2.*t/topt))*               &
        (_ONE_-exp(-self%iv**2*psum**2))
   return
   end function fpz
!EOC

!-----------------------------------------------------------------------

  END MODULE fabm_iow_ergom

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

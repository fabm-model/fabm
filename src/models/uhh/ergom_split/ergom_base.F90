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
! !MODULE: fabm\_uhh\_ergom\_split --- no-phytoplankton version of the biogeochemical model ERGOM
!
! !INTERFACE:
   module fabm_uhh_ergom_split_base
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
   use fabm_uhh_ergom_split_utilities

   implicit none

   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_uhh_ergom_split_base
!
! !REVISION HISTORY:
!  Author(s):
!
! !PUBLIC_DERIVED_TYPES:
   type,extends(type_base_model) :: type_uhh_ergom_split_base
! Variable identifiers
      type (type_state_variable_id)        :: id_de,id_am,id_ni,id_po,id_o2
      type (type_bottom_state_variable_id) :: id_fl
      type (type_dependency_id)            :: id_temp,id_salt
      type (type_horizontal_dependency_id) :: id_wind,id_taub
      type (type_diagnostic_variable_id)   :: id_GPP,id_NCP,id_PPR,id_NPR
! Model parameters
      real(rk) :: sfl_po,sfl_am,sfl_ni,kc,w_de, fl_initial
      real(rk) :: tf,lan,anmx,oan,beta_an,lda,tda,beta_da
      real(rk) :: pvel,sr,s1,s2,s3,s4,a0,a1,a2,lds,lsd,tau_crit,lsa,bsa,ph1,ph2
      real(rk) :: minimum_thsum,depo,tbg
      logical  :: fluff

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
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
   class (type_uhh_ergom_split_base), intent(inout), target :: self
   integer,                intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk)           :: sfl_po=0.0015_rk
   real(rk)           :: sfl_am=0.07_rk
   real(rk)           :: sfl_ni=0.09_rk
   logical            :: fluff=.false.
   real(rk)           :: fl_initial=0.0_rk
   real(rk)           :: w_de=-3.0_rk
   real(rk)           :: kc=0.03_rk
   real(rk)           :: lan=0.1_rk
   real(rk)           :: anmx=0.05_rk
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
   real(rk)           :: depo=0.0_rk
   real(rk)           :: minimum_thsum=0.0_rk
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk
   
   namelist /uhh_ergom_split_base/ &
                        sfl_po,sfl_am,sfl_ni,fluff,    &
                        w_de,kc,depo,   &
                        lan,oan,beta_an,anmx,  &
                        lda,tda,beta_da,pvel,sr,s1,s2,s3,s4,a0,a1,a2,  &
                        lds,lsd,tau_crit,lsa,bsa,ph1,ph2,minimum_thsum
!EOP
!-----------------------------------------------------------------------
!BOC
   
   !  Read the namelist
   if (configunit>0) read(configunit,nml=uhh_ergom_split_base,err=99,end=100)

!  Store parameter values in our own derived type
!  NB! All rates must be provided in values per day,
!  and are converted here to values per second
   call self%get_parameter(self%fl_initial, 'fl_initial', default=fl_initial)
   call self%get_parameter(self%sfl_po, 'sfl_po', default=sfl_po, scale_factor=one_pr_day)
   call self%get_parameter(self%sfl_ni, 'sfl_ni', default=sfl_ni, scale_factor=one_pr_day)
   call self%get_parameter(self%sfl_am, 'sfl_am', default=sfl_am, scale_factor=one_pr_day)
   call self%get_parameter(self%fluff, 'fluff', default=fluff) 
   call self%get_parameter(self%kc, 'kc', default=kc)
   call self%get_parameter(self%lan, 'lan', default=lan, scale_factor=one_pr_day)
   call self%get_parameter(self%anmx, 'anmx', default=anmx, scale_factor=one_pr_day)
   call self%get_parameter(self%oan, 'oan', default=oan)
   call self%get_parameter(self%beta_an, 'beta_an', default=beta_an)
   call self%get_parameter(self%lda, 'lda', default=lda, scale_factor=one_pr_day)
   call self%get_parameter(self%tda, 'tda', default=tda)
   call self%get_parameter(self%beta_da, 'beta_da', default=beta_da)
   call self%get_parameter(self%pvel, 'pvel', default=pvel, scale_factor=one_pr_day)
   call self%get_parameter(self%sr, 'sr', default=sr)
   call self%get_parameter(self%s1, 's1', default=s1)
   call self%get_parameter(self%s2, 's2', default=s2)
   call self%get_parameter(self%s3, 's3', default=s3)
   call self%get_parameter(self%s4, 's4', default=s4)
   call self%get_parameter(self%a0, 'a0', default=a0)
   call self%get_parameter(self%a1, 'a1', default=a1)
   call self%get_parameter(self%a2, 'a2', default=a2)
   call self%get_parameter(self%depo, 'depo', default=depo)
   call self%get_parameter(self%lds, 'lds', default=lds, scale_factor=one_pr_day)
   call self%get_parameter(self%lsd, 'lsd', default=lsd, scale_factor=one_pr_day)
   call self%get_parameter(self%tau_crit, 'tau_crit', default=tau_crit)
   call self%get_parameter(self%lsa, 'lsa', default=lsa, scale_factor=one_pr_day)
   call self%get_parameter(self%bsa, 'bsa', default=bsa)
   call self%get_parameter(self%ph1, 'ph1', default=ph1)
   call self%get_parameter(self%ph2, 'ph2', default=ph2)
   call self%get_parameter(self%minimum_thsum, 'minimum_thsum', default=minimum_thsum)
   call self%get_parameter(self%w_de, 'w_de', default=w_de,scale_factor=one_pr_day)
   
   !  Register state variables

   call self%register_state_variable(self%id_de,'det','mmol n/m**3','detritus',     &
         minimum=0.0_rk,vertical_movement=self%w_de)
   call self%register_state_variable(self%id_am,'amm','mmol n/m**3','ammonium',     &
         minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_ni,'nit','mmol n/m**3','nitrate',      &
         minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_po,'pho','mmol p/m**3','phosphate',    &
         minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_o2,'oxy','mmol o2/m**3','oxygen') 
   if (self%fluff) call self%register_bottom_state_variable(self%id_fl,'flf','mmol n/m**2','fluff', &
          fl_initial,  minimum=0.0_rk)

   ! Register environmental dependencies
! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_wind,standard_variables%wind_speed)
   if (self%fluff) call self%register_dependency(self%id_taub,standard_variables%bottom_stress)
   call self%register_dependency(self%id_taub,standard_variables%bottom_stress) 
   return

 99 call self%fatal_error('uhh_ergom_split_base_init','Error reading namelist uhh_ergom_split_base.')

100 call self%fatal_error('uhh_ergom_split_base_init','Namelist uhh_ergom_split_base was not found.')

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
   class(type_uhh_ergom_split_base), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_

! !LOCAL VARIABLES:
    real(rk)        :: am,ni,po,de,o2,o2pos
    real(rk)        :: temp,llda,llan
    real(rk)        :: wo=30.0_rk,wn=0.1_rk
    real(rk)        :: thopnp,thomnp,thomnm,thsum
    real(rk)        :: o2rhs
!EOP
!-----------------------------------------------------------------------
!BOC
!   Enter spatial_loops (if any)
    _LOOP_BEGIN_

!   Retrieve current (local) state variable values
    _GET_(self%id_de,de) !detritus
    _GET_(self%id_am,am) !ammonium
    _GET_(self%id_ni,ni) !nitrate
    _GET_(self%id_po,po) !phosphate
    _GET_(self%id_o2,o2) !oxygen

!   Retrieve current environmental conditions
    _GET_   (self%id_temp,temp) !local water temperature

    thopnp=th( o2,wo,0.0_rk,1.0_rk)*yy(wn,ni)
    thomnp=th(-o2,wo,0.0_rk,1.0_rk)*yy(wn,ni)
    thomnm=th(-o2,wo,0.0_rk,1.0_rk)*(1.0_rk-yy(wn,ni))
    thsum=thopnp+thomnp+thomnm+self%minimum_thsum
    thopnp=thopnp/thsum
    thomnp=thomnp/thsum
    thomnm=thomnm/thsum

    llda=self%lda*(1.0_rk+self%beta_da*yy(self%tda,temp))
    o2pos = max(o2,0.0_rk)
    llan=th(o2,0.0_rk,0.0_rk,1.0_rk)*o2pos/(self%oan+o2pos)*self%lan*exp(self%beta_an*temp)

   _ADD_SOURCE_(self%id_de,-llda*de)
   _ADD_SOURCE_(self%id_am,llda*de-llan*am-self%anmx*am*thomnp)
   _ADD_SOURCE_(self%id_ni,llan*am-self%s1*llda*de*thomnp)
   _ADD_SOURCE_(self%id_po,self%sr*llda*de)
   o2rhs = -self%s4*(llan*am) - self%s2*(thopnp+thomnm)*llda*de
   
   _ADD_SOURCE_(self%id_o2,o2rhs)

!   Leave spatial loops (if any)
   _LOOP_END_
  
   END subroutine do
!EOC

   subroutine do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   class(type_uhh_ergom_split_base), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_PPDD_
!
!
! !LOCAL VARIABLES:
    real(rk)        :: am,ni,po,de,o2
    real(rk)        :: temp,llda,llan
    real(rk)        :: wo=30.0_rk,wn=0.1_rk
    real(rk)        :: thopnp,thomnp,thomnm,thsum
    real(rk)        :: o2rhs
!EOP
!-----------------------------------------------------------------------
!BOC
!   Enter spatial_loops (if any)
    _LOOP_BEGIN_

!   Retrieve current (local) state variable values
    _GET_(self%id_de,de) !detritus
    _GET_(self%id_am,am) !ammonium
    _GET_(self%id_ni,ni) !nitrate
    _GET_(self%id_po,po) !phosphate
    _GET_(self%id_o2,o2) !oxygen

!   Retrieve current environmental conditions
    _GET_   (self%id_temp,temp) !local water temperature

    thopnp=th( o2,wo,0.0_rk,1.0_rk)*yy(wn,ni)
    thomnp=th(-o2,wo,0.0_rk,1.0_rk)*yy(wn,ni)
    thomnm=th(-o2,wo,0.0_rk,1.0_rk)*(1.0_rk-yy(wn,ni))
    thsum=thopnp+thomnp+thomnm
    thopnp=thopnp/thsum
    thomnp=thomnp/thsum
    thomnm=thomnm/thsum

    llda=self%lda*(1.0_rk+self%beta_da*yy(self%tda,temp))
    llan=th(o2,0.0_rk,0.0_rk,1.0_rk)*o2/(self%oan+o2)*self%lan*exp(self%beta_an*temp)

   _SET_DD_SYM_(self%id_de,self%id_am,llda*de)
   _SET_DD_SYM_(self%id_am,self%id_ni,llan*am)
   _SET_DD_(self%id_am,self%id_am,self%anmx*am*thomnp)
   _SET_DD_(self%id_ni,self%id_ni,self%s1*llda*de*thomnp)
   _SET_PP_(self%id_po,self%id_po,self%sr*llda*de)
   o2rhs = self%s2*(thopnp+thomnm)*llda*de + self%s4*(llan*am)
   
   _SET_DD_(self%id_o2,self%id_o2,o2rhs)

!   Leave spatial loops (if any)
   _LOOP_END_

   END subroutine do_ppdd


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
   class (type_uhh_ergom_split_base), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                     :: de
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   _GET_(self%id_de,de) ! detritus

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*de)

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
   class (type_uhh_ergom_split_base), intent(in) :: self
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

      _ADD_BOTTOM_SOURCE_(self%id_fl,llds*deb-llsd*fl-llsa*fl-th(oxb,wo,0.0_rk,1.0_rk)*llsa*fl)
      _ADD_BOTTOM_FLUX_(self%id_de,-llds*deb+llsd*fl)
      _ADD_BOTTOM_FLUX_(self%id_am,llsa*fl)
      _ADD_BOTTOM_FLUX_(self%id_ni,-self%s1*thomnp*llsa*fl)
      _ADD_BOTTOM_FLUX_(self%id_po,self%sr*(1.0_rk-self%ph1*th(oxb,wo,0.0_rk,1.0_rk)*yy(self%ph2,oxb))*llsa*fl)
      _ADD_BOTTOM_FLUX_(self%id_o2,-(self%s4+self%s2*(thopnp+thomnm))*llsa*fl)
   else

   ! if running without fluff layer, use burial loss of detritus
      _GET_(self%id_de,deb)
      _ADD_BOTTOM_FLUX_(self%id_de,-self%depo*deb*deb)
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
   class (type_uhh_ergom_split_base),intent(in) :: self
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
   _ADD_SURFACE_FLUX_(self%id_o2,flo2)
   else
      flo2 = self%pvel*(self%a0*(self%a1-self%a2*temp)-o2)
   _ADD_SURFACE_FLUX_(self%id_o2,flo2)
   end if

   _ADD_SURFACE_FLUX_(self%id_ni,self%sfl_ni)!/secs_pr_day)
   _ADD_SURFACE_FLUX_(self%id_am,self%sfl_am)!/secs_pr_day)
   _ADD_SURFACE_FLUX_(self%id_po,self%sfl_po)!/secs_pr_day)

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

  END MODULE fabm_uhh_ergom_split_base

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

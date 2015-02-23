#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_npzd --- NPZD biogeochemical model based upon
! Fennel \& Neumann (1996, German Journal of Hydrography),
! with minor modifications by Burchard et al. (2005, Ocean Dynamics)
! taken from GOTM and adapted for FABM by Jorn Bruggeman
!
! !INTERFACE:
   module gotm_npzd
!
! !DESCRIPTION:
! The NPZD (nutrient-phytoplankton-zooplankton-detritus) model described here
! consists of $I=4$ state variables.
! Nutrient uptake (phytoplankton growth) is limited by light and nutrient
! availability, the latter of which is modelled by means
! of Michaelis-Menten kinetics, see eq.\ (\ref{dnp}).
! The half-saturation nutrient concentration $\alpha$ used in this
! formulation has typically a value between 0.2 and 1.5 mmol N\, m$^{-3}$.
! Zooplankton grazing which is limited by the phytoplankton standing stock
! is modelled by means of an Ivlev formulation, see eq.\ (\ref{dpz}).
! All other processes are based on linear first-order kinematics,
! see eqs.\ (\ref{dpn}) - (\ref{dzd}).
! For all details of the NPZD model implemented here,
! see \cite{Burchardetal2005b}.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_gotm_npzd
!     Variable identifiers
      type (type_state_variable_id)        :: id_n,id_p,id_z,id_d
      type (type_state_variable_id)        :: id_dic
      type (type_dependency_id)            :: id_par
      type (type_horizontal_dependency_id) :: id_I_0
      type (type_diagnostic_variable_id)   :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR

!     Model parameters
      real(rk) :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
      real(rk) :: dic_per_n
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
      procedure :: get_light_extinction
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_gotm_npzd), intent(inout), target :: self
   integer,                intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)          :: n_initial
   real(rk)          :: p_initial
   real(rk)          :: z_initial
   real(rk)          :: d_initial
   real(rk)          :: p0
   real(rk)          :: z0
   real(rk)          :: w_p
   real(rk)          :: w_d
   real(rk)          :: kc
   real(rk)          :: i_min
   real(rk)          :: rmax
   real(rk)          :: gmax
   real(rk)          :: iv
   real(rk)          :: alpha
   real(rk)          :: rpn
   real(rk)          :: rzn
   real(rk)          :: rdn
   real(rk)          :: rpdu
   real(rk)          :: rpdl
   real(rk)          :: rzd
   real(rk)          :: dic_per_n
   character(len=64) :: dic_variable

   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   namelist /gotm_npzd/ n_initial,p_initial,z_initial,d_initial,   &
                        p0,z0,w_p,w_d,kc,i_min,rmax,gmax,iv,alpha,rpn,  &
                        rzn,rdn,rpdu,rpdl,rzd,dic_variable,dic_per_n
!EOP
!-----------------------------------------------------------------------
!BOC
   n_initial = 4.5_rk
   p_initial = 0.0_rk
   z_initial = 0.0_rk
   d_initial = 4.5_rk
   p0        = 0.0225_rk
   z0        = 0.0225_rk
   w_p       = -1.0_rk
   w_d       = -5.0_rk
   kc        = 0.03_rk
   i_min     = 25.0_rk
   rmax      = 1.0_rk
   gmax      = 0.5_rk
   iv        = 1.1_rk
   alpha     = 0.3_rk
   rpn       = 0.01_rk
   rzn       = 0.01_rk
   rdn       = 0.003_rk
   rpdu      = 0.02_rk
   rpdl      = 0.1_rk
   rzd       = 0.02_rk
   dic_per_n = 106.0_rk/16.0_rk  ! Redfield C:N
   dic_variable = ''

   ! Read the namelist
   if (configunit>0) read(configunit,nml=gotm_npzd,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%p0,'p0','mmol m-3','background phytoplankton concentration ',default=p0)
   call self%get_parameter(self%z0,'z0','mmol m-3','background zooplankton concentration',default=z0)
   call self%get_parameter(self%kc,'kc','m2 mmol-1','specific light extinction of phytoplankton and detritus',default=kc)
   call self%get_parameter(self%i_min,'i_min','W m-2','minimum light intensity in euphotic zone',default=i_min)
   call self%get_parameter(self%rmax,'rmax','d-1','maximum specific growth rate of phytoplankton',default=rmax,scale_factor=d_per_s)
   call self%get_parameter(self%gmax,'gmax','d-1','maximum specific grazing rate of zooplankton',default=gmax,scale_factor=d_per_s)
   call self%get_parameter(self%iv,'iv','m3 mmol-1','Ivlev grazing constant',default=iv)
   call self%get_parameter(self%alpha,'alpha','mmol m-3','half-saturation nutrient concentration for phytoplankton',default=alpha)
   call self%get_parameter(self%rpn,'rpn','d-1','loss rate of phytoplankton to nutrients',default=rpn,scale_factor=d_per_s)
   call self%get_parameter(self%rzn,'rzn','d-1','loss rate of zooplankton to nutrients',default=rzn,scale_factor=d_per_s)
   call self%get_parameter(self%rdn,'rdn','d-1','detritus remineralization rate',default=rdn,scale_factor=d_per_s)
   call self%get_parameter(self%rpdu,'rpdu','d-1','phytoplankton mortality in euphotic zone',default=rpdu,scale_factor=d_per_s)
   call self%get_parameter(self%rpdl,'rpdl','d-1','phytoplankton mortality below euphotic zone',default=rpdl,scale_factor=d_per_s)
   call self%get_parameter(self%rzd,'rzd','d-1','zooplankton mortality',default=rzd,scale_factor=d_per_s)
   call self%get_parameter(self%dic_per_n,'dic_per_n','-','C:N ratio of biomass',default=dic_per_n)
   call self%get_parameter(w_p,'w_p','m d-1','vertical velocity of phytoplankton (<0 for sinking)',default=w_p, scale_factor=d_per_s)
   call self%get_parameter(w_d,'w_d','m d-1','vertical velocity of detritus  (<0 for sinking)',default=w_d,scale_factor=d_per_s)

   ! Register state variables
   call self%register_state_variable(self%id_n,'nut','mmol m-3','nutrients',    n_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_p,'phy','mmol m-3','phytoplankton',p_initial,minimum=0.0_rk,vertical_movement=w_p)
   call self%register_state_variable(self%id_z,'zoo','mmol m-3','zooplankton',  z_initial,minimum=0.0_rk)
   call self%register_state_variable(self%id_d,'det','mmol m-3','detritus',     d_initial,minimum=0.0_rk,vertical_movement=w_d)

   ! Register the contribution of all state variables to total nitrogen
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_p)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_z)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_d)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   call self%register_state_dependency(self%id_dic,'dic','mmol m-3','total dissolved inorganic carbon',required=.false.)
   if (dic_variable/='') call self%request_coupling(self%id_dic,dic_variable)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_GPP,'GPP','mmol m-3',  'gross primary production',           &
                     output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_NCP,'NCP','mmol m-3',  'net community production',           &
                     output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_PPR,'PPR','mmol m-3 d-1','gross primary production rate',      &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_NPR,'NPR','mmol m-3 d-1','net community production rate',      &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W m-2',    'photosynthetically active radiation',&
                     output=output_time_step_averaged)

   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)

   return

99 call self%fatal_error('gotm_npzd_init','Error reading namelist gotm_npzd.')

100 call self%fatal_error('gotm_npzd_init','Namelist gotm_npzd was not found.')

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
! !DESCRIPTION:
! Seven processes expressed as sink terms are included in this
! conservative model, see eqs.\ (\ref{dnp}) - (\ref{dzd}). \\
!
! Nutrient uptake by phytoplankton:
! \begin{equation}\label{dnp}
! d_{np} = r_{\max}\frac{I_{PAR}}{I_{opt}}
! \exp\left(1-\frac{I_{PAR}}{I_{opt}}\right)
! \frac{c_n}{\alpha+c_n}c_p
! \end{equation}
!
! with
!
! \begin{equation}
! I_{opt}=\max\left(\frac14I_{PAR},I_{\min}\right).
! \end{equation}
!
! Grazing of zooplankton on phytoplankton:
! \begin{equation}\label{dpz}
! d_{pz}=g_{\max}\left(1-\exp\left(-I_v^2c_p^2\right)\right)c_z
! \end{equation}
!
! Phytoplankton excretion:
! \begin{equation}\label{dpn}
! d_{pn} = r_{pn} c_p
! \end{equation}
!
! Zooplankton excretion:
! \begin{equation}\label{dzn}
! d_{zn} = r_{zn} c_z
! \end{equation}
!
! Remineralisation of detritus into nutrients:
! \begin{equation}\label{ddn}
! d_{dn} = r_{dn} c_d
! \end{equation}
!
! Phytoplankton mortality:
! \begin{equation}\label{dpd}
! d_{pd} = r_{pd} c_p
! \end{equation}
!
! Zooplankton mortality:
! \begin{equation}\label{dzd}
! d_{zd} = r_{zd} c_z
! \end{equation}
!
! !INPUT PARAMETERS:
   class (type_gotm_npzd),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                   :: n,p,z,d,par,I_0
   real(rk)                   :: iopt,rpd,primprod,dn
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_n,n) ! nutrient
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_z,z) ! zooplankton
   _GET_(self%id_d,d) ! detritus

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation

   ! Light acclimation formulation based on surface light intensity.
   iopt = max(0.25*I_0,self%I_min)

   ! Loss rate of phytoplankton to detritus depends on local light intensity.
   if (par .ge. self%I_min) then
      rpd = self%rpdu
   else
      rpd = self%rpdl
   end if

   ! Define some intermediate quantities that will be reused multiple times.
   primprod = fnp(self,n,p,par,iopt)
   dn = - primprod + self%rpn*p + self%rzn*z + self%rdn*d

   ! Set temporal derivatives
   _SET_ODE_(self%id_n,dn)
   _SET_ODE_(self%id_p,primprod - fpz(self,p,z) - self%rpn*p - rpd*p)
   _SET_ODE_(self%id_z,fpz(self,p,z) - self%rzn*z - self%rzd*z)
   _SET_ODE_(self%id_d,rpd*p + self%rzd*z - self%rdn*d)

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic,self%dic_per_n*dn)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP ,primprod)
   _SET_DIAGNOSTIC_(self%id_NCP ,primprod - self%rpn*p)
   _SET_DIAGNOSTIC_(self%id_PPR ,primprod*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR ,(primprod - self%rpn*p)*secs_pr_day)

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
   class (type_gotm_npzd), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                     :: p,d
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_d,d) ! detritus

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(self%p0+p+d))

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!
! !DESCRIPTION:
! Seven processes expressed as sink terms are included in this
! conservative model, see eqs.\ (\ref{dnp}) - (\ref{dzd}). \\
!
! Nutrient uptake by phytoplankton:
! \begin{equation}\label{dnp}
! d_{np} = r_{\max}\frac{I_{PAR}}{I_{opt}}
! \exp\left(1-\frac{I_{PAR}}{I_{opt}}\right)
! \frac{c_n}{\alpha+c_n}c_p
! \end{equation}
!
! with
!
! \begin{equation}
! I_{opt}=\max\left(\frac14I_{PAR},I_{\min}\right).
! \end{equation}
!
! Grazing of zooplankton on phytoplankton:
! \begin{equation}\label{dpz}
! d_{pz}=g_{\max}\left(1-\exp\left(-I_v^2c_p^2\right)\right)c_z
! \end{equation}
!
! Phytoplankton excretion:
! \begin{equation}\label{dpn}
! d_{pn} = r_{pn} c_p
! \end{equation}
!
! Zooplankton excretion:
! \begin{equation}\label{dzn}
! d_{zn} = r_{zn} c_z
! \end{equation}
!
! Remineralisation of detritus into nutrients:
! \begin{equation}\label{ddn}
! d_{dn} = r_{dn} c_d
! \end{equation}
!
! Phytoplankton mortality:
! \begin{equation}\label{dpd}
! d_{pd} = r_{pd} c_p
! \end{equation}
!
! Zooplankton mortality:
! \begin{equation}\label{dzd}
! d_{zd} = r_{zd} c_z
! \end{equation}
!
! !INPUT PARAMETERS:
   class (type_gotm_npzd),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                   :: n,p,z,d,par,I_0
   real(rk)                   :: iopt,rpd,dn,primprod
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_n,n) ! nutrient
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_z,z) ! zooplankton
   _GET_(self%id_d,d) ! detritus

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)  ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation

   ! Light acclimation formulation based on surface light intensity.
   iopt = max(0.25*I_0,self%I_min)

   ! Loss rate of phytoplankton to detritus depends on local light intensity.
   if (par .ge. self%I_min) then
      rpd = self%rpdu
   else
      rpd = self%rpdl
   end if

   ! Rate of primary production will be reused multiple times - calculate it once.
   primprod = fnp(self,n,p,par,iopt)

   ! Assign destruction rates to different elements of the destruction matrix.
   ! By assigning with _SET_DD_SYM_(i,j,val) as opposed to _SET_DD_(i,j,val),
   ! assignments to dd(i,j) are automatically assigned to pp(j,i) as well.
   _SET_DD_SYM_(self%id_n,self%id_p,primprod)               ! snp
   _SET_DD_SYM_(self%id_p,self%id_z,fpz(self,p,z))          ! spz
   _SET_DD_SYM_(self%id_p,self%id_n,self%rpn*p)             ! spn
   _SET_DD_SYM_(self%id_z,self%id_n,self%rzn*z)             ! szn
   _SET_DD_SYM_(self%id_d,self%id_n,self%rdn*d)             ! sdn
   _SET_DD_SYM_(self%id_p,self%id_d,rpd*p)                  ! spd
   _SET_DD_SYM_(self%id_z,self%id_d,self%rzd*z)             ! szd

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   dn = - fnp(self,n,p,par,iopt) + self%rpn*p + self%rzn*z + self%rdn*d
   if (_AVAILABLE_(self%id_dic)) _SET_PP_(self%id_dic,self%id_dic,self%dic_per_n*dn)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP,primprod)
   _SET_DIAGNOSTIC_(self%id_NCP,primprod-self%rpn*p)
   _SET_DIAGNOSTIC_(self%id_PPR,primprod*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR,(primprod-self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:
   pure real(rk) function fnp(self,n,p,par,iopt)
!
! !DESCRIPTION:
! Here, the classical Michaelis-Menten formulation for nutrient uptake
! is formulated.
!
! !INPUT PARAMETERS:
   type (type_gotm_npzd), intent(in) :: self
   real(rk), intent(in)         :: n,p,par,iopt
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fnp = self%rmax*par/iopt*exp(1.0_rk-par/iopt)*n/(self%alpha+n)*(p+self%p0)

   end function fnp
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure real(rk) function fpz(self,p,z)
!
! !DESCRIPTION:
! Here, the classical Ivlev formulation for zooplankton grazing on
! phytoplankton is formulated.
!
! !INPUT PARAMETERS:
   type (type_gotm_npzd), intent(in) :: self
   real(rk), intent(in)         :: p,z
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fpz = self%gmax*(1.0_rk-exp(-self%iv*self%iv*p*p))*(z+self%z0)

   end function fpz
!EOC

!-----------------------------------------------------------------------

   end module gotm_npzd

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

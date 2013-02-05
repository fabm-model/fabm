#include "fabm_driver.h"

#ifdef _FABM_F2003_

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_npzd --- Fortran 2003 version of NPZD biogeochemical model based upon
! Fennel \& Neumann (1996, German Journal of Hydrography),
! with minor modifications by Burchard et al. (2005, Ocean Dynamics)
! taken from GOTM and adapted for FABM by Jorn Bruggeman
!
! !INTERFACE:
   module fabm_examples_npzd_f2003
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
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public examples_npzd_f2003_create, type_examples_npzd_f2003
!
! !PRIVATE DATA MEMBERS:
   REALTYPE, parameter :: secs_pr_day = 86400.
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_examples_npzd_f2003
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_      :: id_n,id_p,id_z,id_d
      _TYPE_STATE_VARIABLE_ID_      :: id_dic
      _TYPE_DEPENDENCY_ID_          :: id_par
      _TYPE_HORIZONTAL_DEPENDENCY_ID_  :: id_I_0
      _TYPE_DIAGNOSTIC_VARIABLE_ID_ :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR
      _TYPE_CONSERVED_QUANTITY_ID_  :: id_totN

!     Model parameters
      REALTYPE :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
      REALTYPE :: dic_per_n
      logical  :: use_dic

      contains

!     Model procedures
      procedure :: do
      procedure :: do_ppdd
      procedure :: get_light_extinction
      procedure :: get_conserved_quantities

   end type type_examples_npzd_f2003
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   function examples_npzd_f2003_create(configunit,name,parent) result(self)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   integer,                          intent(in)    :: configunit
   character(len=*),                 intent(in)    :: name
   class (type_base_model),target,   intent(inout) :: parent
   class (type_examples_npzd_f2003), pointer       :: self
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                  :: n_initial=4.5
   REALTYPE                  :: p_initial=0.
   REALTYPE                  :: z_initial=0.
   REALTYPE                  :: d_initial=4.5
   REALTYPE                  :: p0=0.0225
   REALTYPE                  :: z0=0.0225
   REALTYPE                  :: w_p=-1.157407e-05
   REALTYPE                  :: w_d=-5.787037e-05
   REALTYPE                  :: kc=0.03
   REALTYPE                  :: i_min=25.
   REALTYPE                  :: rmax=1.157407e-05
   REALTYPE                  :: gmax=5.787037e-06
   REALTYPE                  :: iv=1.1
   REALTYPE                  :: alpha=0.3
   REALTYPE                  :: rpn=1.157407e-07
   REALTYPE                  :: rzn=1.157407e-07
   REALTYPE                  :: rdn=3.472222e-08
   REALTYPE                  :: rpdu=2.314814e-07
   REALTYPE                  :: rpdl=1.157407e-06
   REALTYPE                  :: rzd=2.314814e-07
   REALTYPE                  :: dic_per_n=106.d0/16.d0
   character(len=64)         :: dic_variable=''

   namelist /examples_npzd_f2003/ n_initial,p_initial,z_initial,d_initial,   &
                        p0,z0,w_p,w_d,kc,i_min,rmax,gmax,iv,alpha,rpn,  &
                        rzn,rdn,rpdu,rpdl,rzd,dic_variable,dic_per_n
!EOP
!-----------------------------------------------------------------------
!BOC
   allocate(self)
   call self%initialize(name,parent)

   ! Read the namelist
   read(configunit,nml=examples_npzd_f2003,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%p0    = p0
   self%z0    = z0
   self%kc    = kc
   self%i_min = i_min
   self%rmax  = rmax/secs_pr_day
   self%gmax  = gmax/secs_pr_day
   self%iv    = iv
   self%alpha = alpha
   self%rpn  = rpn /secs_pr_day
   self%rzn  = rzn /secs_pr_day
   self%rdn  = rdn /secs_pr_day
   self%rpdu = rpdu/secs_pr_day
   self%rpdl = rpdl/secs_pr_day
   self%rzd  = rzd /secs_pr_day
   self%dic_per_n = dic_per_n

   ! Register state variables
   call self%register_state_variable(self%id_n,'nut','mmol/m**3','nutrients',     &
                                    n_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_p,'phy','mmol/m**3','phytoplankton', &
                                    p_initial,minimum=_ZERO_,vertical_movement=w_p/secs_pr_day)
   call self%register_state_variable(self%id_z,'zoo','mmol/m**3','zooplankton', &
                                    z_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_d,'det','mmol/m**3','detritus', &
                                    d_initial,minimum=_ZERO_,vertical_movement=w_d/secs_pr_day)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   self%use_dic = dic_variable.ne.''
   if (self%use_dic) call self%register_state_dependency(self%id_dic,dic_variable)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_GPP,'GPP','mmol/m**3',  'gross primary production',           &
                     time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_NCP,'NCP','mmol/m**3',  'net community production',           &
                     time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_PPR,'PPR','mmol/m**3/d','gross primary production rate',      &
                     time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_NPR,'NPR','mmol/m**3/d','net community production rate',      &
                     time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2',    'photosynthetically active radiation',&
                     time_treatment=time_treatment_averaged)

   ! Register conserved quantities
   call self%register_conserved_quantity(self%id_totN,'N','mmol/m**3','nitrogen')

   ! Register environmental dependencies
   call self%register_dependency(self%id_par,varname_par)
   call self%register_dependency(self%id_I_0,varname_par_sf)

   return

99 call fatal_error('examples_npzd_f2003_init','Error reading namelist examples_npzd_f2003.')

100 call fatal_error('examples_npzd_f2003_init','Namelist examples_npzd_f2003 was not found.')

   end function examples_npzd_f2003_create
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine do(self,_FABM_ARGS_DO_RHS_)
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
   class (type_examples_npzd_f2003),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: n,p,z,d,par,I_0
   REALTYPE                   :: iopt,rpd,primprod,dn
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

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
   if (self%use_dic) then
      _SET_ODE_(self%id_dic,self%dic_per_n*dn)
   end if

   ! Export diagnostic variables
   _SET_DIAG_(self%id_dPAR,par)
   _SET_DIAG_(self%id_GPP ,primprod)
   _SET_DIAG_(self%id_NCP ,primprod - self%rpn*p)
   _SET_DIAG_(self%id_PPR ,primprod*secs_pr_day)
   _SET_DIAG_(self%id_NPR ,(primprod - self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_f2003), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: p,d
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_d,d) ! detritus

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(self%p0+p+d))

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   subroutine get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_f2003), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: n,p,z,d
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_n,n) ! nutrient
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_z,z) ! zooplankton
   _GET_(self%id_d,d) ! detritus

   ! Total nutrient is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totN,n+p+z+d)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine do_ppdd(self,_FABM_ARGS_DO_PPDD_)
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
   class (type_examples_npzd_f2003),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: n,p,z,d,par,I_0
   REALTYPE                   :: iopt,rpd,dn,primprod
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_n,n) ! nutrient
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_z,z) ! zooplankton
   _GET_(self%id_d,d) ! detritus

   ! Retrieve current environmental conditions.
   _GET_   (self%id_par,par)  ! local photosynthetically active radiation
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
   if (self%use_dic) _SET_PP_(self%id_dic,self%id_dic,self%dic_per_n*dn)

   ! Export diagnostic variables
   _SET_DIAG_(self%id_dPAR,par)
   _SET_DIAG_(self%id_GPP,primprod)
   _SET_DIAG_(self%id_NCP,primprod-self%rpn*p)
   _SET_DIAG_(self%id_PPR,primprod*secs_pr_day)
   _SET_DIAG_(self%id_NPR,(primprod-self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:
   pure REALTYPE function fnp(self,n,p,par,iopt)
!
! !DESCRIPTION:
! Here, the classical Michaelis-Menten formulation for nutrient uptake
! is formulated.
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_f2003), intent(in) :: self
   REALTYPE, intent(in)         :: n,p,par,iopt
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fnp = self%rmax*par/iopt*exp(_ONE_-par/iopt)*n/(self%alpha+n)*(p+self%p0)

   end function fnp
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure REALTYPE function fpz(self,p,z)
!
! !DESCRIPTION:
! Here, the classical Ivlev formulation for zooplankton grazing on
! phytoplankton is formulated.
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_f2003), intent(in) :: self
   REALTYPE, intent(in)         :: p,z
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fpz = self%gmax*(_ONE_-exp(-self%iv*self%iv*p*p))*(z+self%z0)

   end function fpz
!EOC

!-----------------------------------------------------------------------

   end module fabm_examples_npzd_f2003

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#endif

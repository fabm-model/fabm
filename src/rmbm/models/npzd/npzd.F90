!$Id$
#include "rmbm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: rmbm_npzd --- 0D NPZD biogeochemical model
!
! !INTERFACE:
   module rmbm_npzd
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
   use rmbm_types
   use rmbm_driver

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_npzd, npzd_init, npzd_do, npzd_do_ppdd, &
          npzd_get_light_extinction, npzd_get_conserved_quantities
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!
! !PUBLIC DERIVED TYPES:
   type type_npzd
      ! State variable identifiers
      integer  :: id_n,id_p,id_z,id_d

      ! Environmental variable identifiers
      integer  :: id_par,id_I_0
      
      ! Dependencies
      integer  :: id_dic

      ! Diagnostic variable identifiers
      integer  :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR

      ! Conserved quantity identifiers
      integer  :: id_totN
      
      ! Model parameters
      REALTYPE :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
      REALTYPE :: dic_per_n
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine npzd_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio\_npzd.nml} is read and 
!  various variables (rates and settling velocities) 
!  are transformed into SI units.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_npzd),      intent(out)   :: self
   type (type_model_info),intent(inout) :: modelinfo
   integer,               intent(in )   :: namlst
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

   REALTYPE, parameter :: secs_pr_day = 86400.
   namelist /bio_npzd_nml/ n_initial,p_initial,z_initial,d_initial,   &
                           p0,z0,w_p,w_d,kc,i_min,rmax,gmax,iv,alpha,rpn,  &
                           rzn,rdn,rpdu,rpdl,rzd,dic_variable,dic_per_n
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=bio_npzd_nml,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day, and are converted here to values per second.
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
   
   ! State variables
   self%id_n = register_state_variable(modelinfo,'nut','mmol/m**3','nutrients',     &
                                    n_initial,minimum=_ZERO_,no_river_dilution=.true.)
   self%id_p = register_state_variable(modelinfo,'phy','mmol/m**3','phytoplankton', &
                                    p_initial,minimum=_ZERO_,vertical_movement=w_p/secs_pr_day, &
                                    mussels_inhale=.true.)
   self%id_z = register_state_variable(modelinfo,'zoo','mmol/m**3','zooplankton', &
                                    z_initial,minimum=_ZERO_, &
                                    mussels_inhale=.true.)
   self%id_d = register_state_variable(modelinfo,'det','mmol/m**3','detritus', &
                                    d_initial,minimum=_ZERO_,vertical_movement=w_d/secs_pr_day, &
                                    mussels_inhale=.true.)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   self%id_dic = -1
   if (dic_variable.ne.'') then
      self%id_dic = register_state_dependency(modelinfo,dic_variable)
      if (self%id_dic.eq.-1) call fatal_error('init_bio_npzd','Cannot locate external DIC variable '//dic_variable)
   end if

   ! Diagnostic variables
   self%id_GPP  = register_diagnostic_variable(modelinfo,'GPP','mmol/m**3',  'gross primary production',           3)
   self%id_NCP  = register_diagnostic_variable(modelinfo,'NCP','mmol/m**3',  'net community production',           3)
   self%id_PPR  = register_diagnostic_variable(modelinfo,'PPR','mmol/m**3/d','gross primary production rate',      2)
   self%id_NPR  = register_diagnostic_variable(modelinfo,'NPR','mmol/m**3/d','net community production rate',      2)
   self%id_dPAR = register_diagnostic_variable(modelinfo,'PAR','W/m**2',     'photosynthetically active radiation',2)
   
   ! Conserved quantities
   self%id_totN = register_conserved_quantity(modelinfo,'N','mmol/m**3','nitrogen')
   
   ! Environmental dependencies
   self%id_par = register_dependency(modelinfo, varname_par)
   self%id_I_0 = register_dependency(modelinfo, varname_par_sf, shape=shape2d)

   return

99 call fatal_error('init_bio_npzd','I could not read namelist bio_npzd_nml')
   
   end subroutine npzd_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   _PURE function npzd_get_light_extinction(self,environment,LOCATION) result(extinction)
!
! !INPUT PARAMETERS:
   type (type_npzd),       intent(in) :: self
   type (type_environment),intent(in) :: environment
   LOCATION_TYPE,          intent(in) :: LOCATION
   REALTYPE                           :: extinction
   
   REALTYPE                     :: p,d
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Retrieve current (local) state variable values.
   p = environment%state(self%id_p)%data INDEX_LOCATION
   d = environment%state(self%id_d)%data INDEX_LOCATION
   
   ! Self-shading with explicit contribution from background phytoplankton concentration.
   extinction = self%kc*(self%p0+p+d)

   end function npzd_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   _PURE subroutine npzd_get_conserved_quantities(self,environment,LOCATION,sums)
!
! !INPUT PARAMETERS:
   type (type_npzd),       intent(in)    :: self
   type (type_environment),intent(in)    :: environment
   LOCATION_TYPE,          intent(in)    :: LOCATION
   REALTYPE,               intent(inout) :: sums(:)
   
   REALTYPE              :: n,p,z,d
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Retrieve current (local) state variable values.
   n = environment%state(self%id_n)%data INDEX_LOCATION
   p = environment%state(self%id_p)%data INDEX_LOCATION
   z = environment%state(self%id_z)%data INDEX_LOCATION
   d = environment%state(self%id_d)%data INDEX_LOCATION
   
   ! Total nutrient is simply the sum of all variables.
   sums(self%id_totN) = n+p+z+d

   end subroutine npzd_get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:
   _PURE REALTYPE function fnp(self,n,p,par,iopt)
!
! !DESCRIPTION:
! Here, the classical Michaelis-Menten formulation for nutrient uptake
! is formulated.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_npzd), intent(in) :: self
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
   _PURE REALTYPE function fpz(self,p,z)
!
! !DESCRIPTION:
! Here, the classical Ivlev formulation for zooplankton grazing on 
! phytoplankton is formulated.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_npzd), intent(in) :: self
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
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   _PURE subroutine npzd_do(self,environment,LOCATION,rhs,diag)
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
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_npzd),       intent(in) :: self
   type (type_environment),intent(in) :: environment
   LOCATION_TYPE,          intent(in) :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, dimension(:), intent(inout) :: rhs,diag
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: n,p,z,d,par,I_0,dn
   REALTYPE                   :: iopt
   REALTYPE                   :: rpd
   REALTYPE, parameter :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Retrieve current (local) state variable values.
   n = environment%state(self%id_n)%data INDEX_LOCATION
   p = environment%state(self%id_p)%data INDEX_LOCATION
   z = environment%state(self%id_z)%data INDEX_LOCATION
   d = environment%state(self%id_d)%data INDEX_LOCATION
   
   ! Retrieve current (local) environmental conditions.
   par = environment%var3d(self%id_par)%data INDEX_LOCATION
   I_0 = environment%var2d(self%id_I_0)%data INDEX_LOCATION_HZ
   
   ! Light acclimation formulation based on surface light intensity.
   iopt = max(0.25*I_0,self%I_min)

   if (par .ge. self%I_min) then
      rpd = self%rpdu
   else
      rpd = self%rpdl
   end if
   
   ! Note: the temporal derivatives are incremented or decremented, rahter than set.
   ! This is IMPORTANT: other biogeochemical models might already have provided additional
   ! sink and source terms for these variables in the rhs arrays.
   dn = - fnp(self,n,p,par,iopt) + self%rpn*p + self%rzn*z + self%rdn*d
   rhs(self%id_n) = rhs(self%id_n) + dn
   rhs(self%id_p) = rhs(self%id_p) + fnp(self,n,p,par,iopt) - fpz(self,p,z) - self%rpn*p - rpd*p
   rhs(self%id_z) = rhs(self%id_z) + fpz(self,p,z) - self%rzn*z - self%rzd*z
   rhs(self%id_d) = rhs(self%id_d) + rpd*p + self%rzd*z - self%rdn*d

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   if (self%id_dic.ne.-1) rhs(self%id_dic) = rhs(self%id_dic) + self%dic_per_n*dn

   ! Export diagnostic variables
   if (self%id_dPAR.ne.id_not_used) diag(self%id_dPAR) = par
   if (self%id_GPP .ne.id_not_used) diag(self%id_GPP)  = fnp(self,n,p,par,iopt)
   if (self%id_NCP .ne.id_not_used) diag(self%id_NCP)  = fnp(self,n,p,par,iopt) - self%rpn*p
   if (self%id_PPR .ne.id_not_used) diag(self%id_PPR)  = diag(self%id_GPP)*secs_pr_day
   if (self%id_NPR .ne.id_not_used) diag(self%id_NPR)  = diag(self%id_NCP)*secs_pr_day

   end subroutine npzd_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model exporting production/destruction matrices
!
! !INTERFACE:
   _PURE subroutine npzd_do_ppdd(self,environment,LOCATION,pp,dd,diag)
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
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_npzd),       intent(in) :: self
   type (type_environment),intent(in) :: environment
   LOCATION_TYPE,          intent(in) :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout),dimension(:,:) :: pp,dd
   REALTYPE, intent(inout),dimension(:)   :: diag
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: n,p,z,d,par,I_0,dn
   REALTYPE                   :: iopt
   REALTYPE                   :: rpd
   REALTYPE, parameter :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Retrieve current (local) state variable values.
   n = environment%state(self%id_n)%data INDEX_LOCATION
   p = environment%state(self%id_p)%data INDEX_LOCATION
   z = environment%state(self%id_z)%data INDEX_LOCATION
   d = environment%state(self%id_d)%data INDEX_LOCATION
   
   ! Retrieve current (local) environmental conditions.
   par = environment%var3d(self%id_par)%data INDEX_LOCATION
   I_0 = environment%var2d(self%id_I_0)%data INDEX_LOCATION_HZ
   
   ! Light acclimation formulation based on surface light intensity.
   iopt = max(0.25*I_0,self%I_min)

   if (par .ge. self%I_min) then
      rpd = self%rpdu
   else
      rpd = self%rpdl
   end if
   
   dd(self%id_n,self%id_p) = fnp(self,n,p,par,iopt)  ! snp
   dd(self%id_p,self%id_z) = fpz(self,p,z)           ! spz
   dd(self%id_p,self%id_n) = self%rpn*p              ! spn
   dd(self%id_z,self%id_n) = self%rzn*z              ! szn
   dd(self%id_d,self%id_n) = self%rdn*d              ! sdn
   dd(self%id_p,self%id_d) =      rpd*p              ! spd
   dd(self%id_z,self%id_d) = self%rzd*z              ! szd
   
   ! Mirror destruction rates in production rate matrix (each conversion is fully conservative)
   pp(self%id_p,self%id_n) = dd(self%id_n,self%id_p)
   pp(self%id_z,self%id_p) = dd(self%id_p,self%id_z)
   pp(self%id_n,self%id_p) = dd(self%id_p,self%id_n)
   pp(self%id_n,self%id_z) = dd(self%id_z,self%id_n)
   pp(self%id_n,self%id_d) = dd(self%id_d,self%id_n)
   pp(self%id_d,self%id_p) = dd(self%id_p,self%id_d)
   pp(self%id_d,self%id_z) = dd(self%id_z,self%id_d)

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   dn = - fnp(self,n,p,par,iopt) + self%rpn*p + self%rzn*z + self%rdn*d
   if (self%id_dic.ne.-1) pp(self%id_dic,self%id_dic) = pp(self%id_dic,self%id_dic) + self%dic_per_n*dn

   ! Export diagnostic variables
   if (self%id_dPAR.ne.id_not_used) diag(self%id_dPAR) = par
   if (self%id_GPP .ne.id_not_used) diag(self%id_GPP)  = dd(self%id_n,self%id_p)
   if (self%id_NCP .ne.id_not_used) diag(self%id_NCP)  = dd(self%id_n,self%id_p)-pp(self%id_n,self%id_p)
   if (self%id_PPR .ne.id_not_used) diag(self%id_PPR)  = diag(self%id_GPP)*secs_pr_day
   if (self%id_NPR .ne.id_not_used) diag(self%id_NPR)  = diag(self%id_NCP)*secs_pr_day

   end subroutine npzd_do_ppdd
!EOC


!-----------------------------------------------------------------------

   end module rmbm_npzd

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

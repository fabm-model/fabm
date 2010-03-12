!$Id$
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Heat and momemtum fluxes according to Kondo \label{sec:kondo}
!
! !INTERFACE:
   subroutine kondo(sst,airt,u10,v10,precip,evap,taux,tauy,qe,qh)
!
! !DESCRIPTION:
!  Based on the model sea surface temperature, the wind vector
!  at 10 m height, the air pressure at 2 m, the dry air
!  temperature and the air pressure at 2 m, and the relative
!  humidity (either directly given or recalculated from the
!  wet bulb or the dew point temperature),
!  this routine first computes the transfer coefficients for the surface
!  momentum flux vector, $(\tau_x^s,\tau_y^s)$ ($c_{dd}$), 
!  the latent heat flux, $Q_e$, ($c_{ed}$)
!  and the sensible heat flux, $Q_h$, 
!  ($c_{hd}$) heat flux according to the \cite{Kondo75}
!  bulk formulae. Afterwards, these fluxes are calculated according 
!  to the following formulae:
!
!  \begin{equation}
!  \begin{array}{rcl}
!  \tau_x^s &=& c_{dd} \rho_a W_x W \\ \\
!  \tau_y^s &=& c_{dd} \rho_a W_y W \\ \\
!  Q_e &=& c_{ed} L \rho_a W (q_s-q_a) \\ \\
!  Q_h &=& c_{hd} C_{pa} \rho_a W (T_w-T_a)
!  \end{array}
!  \end{equation}
!
!  with the air density $\rho_a$, the wind speed at 10 m, $W$,
!  the $x$- and the $y$-component of the wind velocity vector,
!  $W_x$ and $W_y$, respectively, the specific evaporation heat of sea water,
!  $L$, the specific saturation humidity, $q_s$, the actual
!  specific humidity $q_a$, the specific heat capacity of air at constant
!  pressure, $C_{pa}$, the sea surface temperature, $T_w$ and the
!  dry air temperature, $T_a$.
!
! !USES:
   use airsea_variables, only: kelvin,const06,rgas,rho_0
   use airsea_variables, only: qs,qa,rhoa
   use airsea_variables, only: cpa,cpw
   use airsea, only: rain_impact,calc_evaporation
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: sst,airt,u10,v10,precip
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)             :: evap
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: taux,tauy,qe,qh
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard and Karsten Bolding
!
!  $Log: kondo.F90,v $
!  Revision 1.2  2007-12-21 12:38:03  kb
!  added precip/evap to kondo + cleaned
!
!  Revision 1.1  2007-09-25 10:06:10  kbk
!  modularized the airsea module - added Fairall method
!
! !LOCAL VARIABLES:
   REALTYPE                  :: w,L
   REALTYPE                  :: s,s0
   REALTYPE                  :: ae_d,be_d,pe_d
   REALTYPE                  :: ae_h,be_h,ce_h,pe_h
   REALTYPE                  :: ae_e,be_e,ce_e,pe_e
   REALTYPE                  :: x,x1,x2,x3
   REALTYPE                  :: ta,ta_k,tw,tw_k
   REALTYPE                  :: cdd,chd,ced
   REALTYPE                  :: tmp,rainfall,cd_rain
   REALTYPE, parameter       :: eps=1.0e-12
!EOP
!-----------------------------------------------------------------------
!BOC
   w = sqrt(u10*u10+v10*v10)
   L = (2.5-0.00234*sst)*1.e6

   if (sst .lt. 100.) then
      tw  = sst
      tw_k= sst+kelvin
   else
      tw  = sst-kelvin
      tw_k= sst
   end if

   if (airt .lt. 100.) then
      ta_k  = airt + kelvin
      ta = airt
   else
      ta  = airt - kelvin
      ta_k = airt
   end if

!  Stability
   s0=0.25*(sst-airt)/(w+1.0e-10)**2
   s=s0*abs(s0)/(abs(s0)+0.01)

!  Transfer coefficient for heat and momentum
   if (w .lt. 2.2) then
      ae_d=0.0;   be_d=1.08;                  pe_d=-0.15;
      ae_h=0.0;   be_h=1.185;  ce_h=0.0;      pe_h=-0.157;
      ae_e=0.0;   be_e=1.23;   ce_e=0.0;      pe_e=-0.16;
   else if (w .lt. 5.0) then
      ae_d=0.771; be_d=0.0858;                pe_d=1.0;
      ae_h=0.927; be_h=0.0546; ce_h=0.0;      pe_h=1.0;
      ae_e=0.969; be_e=0.0521; ce_e=0.0;      pe_e=1.0;
   else if (w .lt. 8.0) then
      ae_d=0.867; be_d=0.0667;                pe_d=1.0;
      ae_h=1.15;  be_h=0.01;   ce_h=0.0;      pe_h=1.0;
      ae_e=1.18;  be_e=0.01;   ce_e=0.0;      pe_e=1.0;
   else if (w .lt. 25.0) then
      ae_d=1.2;   be_d=0.025;                 pe_d=1.0;
      ae_h=1.17;  be_h=0.0075; ce_h=-0.00045; pe_h=1.0;
      ae_e=1.196; be_e=0.008;  ce_e=-0.0004;  pe_e=1.0
   else
      ae_d=0.0;   be_d=0.073;                 pe_d=1.0;
      ae_h=1.652; be_h=-0.017; ce_h=0.0;      pe_h=1.0;
      ae_e=1.68;  be_e=-0.016; ce_e=0;        pe_e=1.0;
   end if

   cdd=(ae_d+be_d*exp(pe_d*log(w+eps)))*1.0e-3
   chd=(ae_h+be_h*exp(pe_h*log(w+eps))+ce_h*(w-8.0)**2)*1.0e-3
   ced=(ae_e+be_e*exp(pe_e*log(w+eps))+ce_e*(w-8.0)**2)*1.0e-3

   if(s .lt. 0.) then
      if (s .gt. -3.3) then
         x = 0.1+0.03*s+0.9*exp(4.8*s)
      else
         x = 0.0
      end if
      cdd=x*cdd
      chd=x*chd
      ced=x*ced
   else
      cdd=cdd*(1.0+0.47*sqrt(s))
      chd=chd*(1.0+0.63*sqrt(s))
      ced=ced*(1.0+0.63*sqrt(s))
   end if

   qh=-ced*L*rhoa*w*(qs-qa)            ! latent

   qe=-chd*cpa*rhoa*w*(sst-airt)       ! sensible

!  compute sensible heatflux correction due to rain fall
   if (rain_impact) then
!     units of qs and qa - should be kg/kg
      rainfall=precip * 1000. ! (convert from m/s to kg/m2/s)
      x1 = 2.11e-5*(ta_k/kelvin)**1.94
      x2 = 0.02411*(1.0+ta*(3.309e-3-1.44e-6*ta))/(rhoa*cpa)
      x3 = qa * L /(rgas * ta_K * ta_K)
      cd_rain = 1.0/(1.0+const06*(x3*L*x1)/(cpa*x2))
      cd_rain = cd_rain*cpw*((tw-ta) + (qs-qa)*L/cpa)
      qe = qe - rainfall * cd_rain
   end if

!  calculation of evaporation/condensation in m/s
   if (rain_impact .and. calc_evaporation) then
!     ced from latent heatflux for moisture flux
      evap = rhoa/rho_0*ced*w*(qa-qs)
   else
      evap = _ZERO_
   end if

   tmp = cdd*rhoa*w
   taux = tmp*u10
   tauy = tmp*v10

!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
   if ( rain_impact ) then
      tmp  = 0.85 * rainfall * w
      x=u10
      taux  = taux + tmp * sign(_ONE_,x)
      x=v10
      tauy  = tauy + tmp * sign(_ONE_,x)
   end if

   return
   end subroutine kondo
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

!$Id$
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Heat and momentum fluxes according to Fairall et al.
!
! !INTERFACE:
   subroutine fairall(sst,airt,u10,v10,precip,evap,taux,tauy,qe,qh)
!
! !DESCRIPTION:
!  The surface momentum flux vector, $(\tau_x^s,\tau_y^s)$,
!  in [N\,m$^{-2}$],
!  the latent heat flux, $Q_e$, 
!  and the sensible heat flux, $Q_h$, both in [W\,m$^{-2}$]
!  are calculated here according to the \cite{Fairalletal96a} bulk
!  formulae, which are build on the Liu-Katsaros-Businger
!  (\cite{Liuetal79}) method. 
!  Cool skin and warm layer effects are considered according to the
!  suggestions of \cite{Fairalletal96b}.
!
!  The air temperature {\tt airt} and the sea surface temperature
!  {\tt sst} may be given in Kelvin or Celcius: 
!  if they are $>$ 100 - Kelvin is assumed. 
!  
!  This piece of code has been adapted from the COARE code originally 
!  written by David Rutgers and Frank Bradley - see
!  http://www.coaps.fsu.edu/COARE/flux\_algor/flux.html.
!
! !USES:
   use airsea_variables, only: kelvin,const06,rgas,rho_0,g,rho_0,kappa
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
!  Original author(s): Adolf Stips
!
!  $Log: fairall.F90,v $
!  Revision 1.8  2009-01-07 07:25:38  kb
!  fixed various compilation warnings found by gfortran
!
!  Revision 1.7  2008-07-07 09:05:08  lars
!  corrected typo in documentation
!
!  Revision 1.6  2008-04-08 16:09:00  kb
!  assure valid qh and qe under all circumstance - Bruggeman, Stips
!
!  Revision 1.5  2008-01-02 15:30:44  kb
!  added link to Fairall page
!
!  Revision 1.4  2007-12-21 12:38:03  kb
!  added precip/evap to kondo + cleaned
!
!  Revision 1.3  2007-12-19 10:41:20  kb
!  fixed m/s --> kg/m2/s conversion bug - Stips
!
!  Revision 1.2  2007-10-02 10:14:08  kbk
!  fixed rhoa calculation - rgas in airsea_variables module
!
!  Revision 1.1  2007-09-25 10:06:10  kbk
!  modularized the airsea module - added Fairall method
!

! !DEFINED PARAMETERS:
!  Fairall LKB roughness Reynolds number to Von Karman
   REALTYPE,parameter        :: fdg = 1.0          ! non-dimensional

!  Beta parameter evaluated from Fairall low windspeed turbulence data.
   REALTYPE,parameter        :: beta = 1.2         ! non-dimensional

!  Zabl      Height (m) of atmospheric boundary layer.
   REALTYPE,parameter        :: Zabl = 600.0       ! in [m]

   REALTYPE, parameter       :: r3 = 1.0/3.0
!
!  Liu et al. (1979) look-up table coefficients to compute roughness
!  Reynolds number for temperature (rt) and moisture (rq) as a
!  function of wind Reynolds number (rr):
!
!       rt = Liu_a(:,1) * Rr   ** Liu_b(:,1)    temperature
!       rq = Liu_a(:,2) * Rr   ** Liu_b(:,2)    moisture
!
   REALTYPE,parameter, dimension(8,2) :: Liu_a = reshape ( &
                 (/ 0.177,  1.376,    1.026,      1.625,   &
                    4.661, 34.904, 1667.190, 588000.0,     &
                    0.292,  1.808,    1.393,      1.956,   &
                    4.994, 30.709, 1448.680, 298000.0 /),  &
                 (/ 8, 2 /) )

   REALTYPE,parameter, dimension(8,2) :: Liu_b = reshape ( &
                 (/  0.0,    0.929, -0.599, -1.018,        &
                    -1.475, -2.067, -2.907, -3.935,        &
                     0.0,    0.826, -0.528, -0.870,        &
                    -1.297, -1.845, -2.682, -3.616 /),     &
                 (/ 8, 2 /) )

   REALTYPE,parameter, dimension(9) :: Liu_Rr =            &
                 (/    0.0,  0.11,   0.825,   3.0,         &
                      10.0, 30.0,  100.0,   300.0,         &
                    1000.0 /)
!
!  Height (m) of surface air temperature measurement.
   REALTYPE, parameter       ::  zt= 2.0
!  Height (m) of surface air humidity measurement
   REALTYPE, parameter       ::  zq= 2.0
!  Height (m) of surface winds measurement
   REALTYPE, parameter       ::  zw= 10.0
   integer,  parameter       :: itermax = 20
#ifdef GUSTINESS
   REALTYPE, parameter       :: wgust=0.2
#else
   REALTYPE, parameter       :: wgust=0.0
#endif

   REALTYPE,external         :: psi

! !LOCAL VARIABLES:
   REALTYPE                  :: tmp,cff,wgus
   REALTYPE                  :: L
   REALTYPE                  :: Cd
   REALTYPE                  :: ta,ta_k,tw,tw_k
   integer                   :: ier,iter,k
   REALTYPE                  :: vis_air
   REALTYPE                  :: tpsi,qpsi,wpsi,ZWoL,oL,ZToL,ZQoL,ZoW,ZoT, ZoQ
   REALTYPE                  :: Wstar,Tstar, Qstar, delQ, delT, rr,rt,rq
   REALTYPE                  :: TVstar,Bf, upvel,delw,Wspeed, w
   REALTYPE                  :: ri,cd_rain
   REALTYPE                  :: x1,x2,x3
   REALTYPE                  :: x
   REALTYPE                  :: rainfall
   REALTYPE, parameter       :: eps=1.0e-12
!EOP
!-----------------------------------------------------------------------
!BOC
   evap = _ZERO_
   w = sqrt(u10*u10+v10*v10)

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

!
!  Initialize.
!
   qe=_ZERO_
   qh=_ZERO_
   delw=sqrt(w*w+wgust*wgust)
   if (delw .ne. 0.0) then
!-----------------------------------------------------------------------
!     Compute Monin-Obukhov similarity parameters for wind (Wstar),
!     heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------

!     Kinematic viscosity of dry air (m2/s), Andreas (1989).
      vis_air=1.326e-5*(1.0+ta*(6.542e-3+ta*(8.301e-6-4.84e-9*ta)))

!     Compute latent heat of vaporization (J/kg) at sea surface
      L = (2.501-0.00237*tw)*1.e6
!
!     Assume that wind is measured relative to sea surface and include
!     gustiness.
!     Initialize.
      ier = 0.0
      delq=qa-qs
      delt=ta-tw

!     Initial guesses for Monin-Obukhov similarity scales.
      ZWoL=0.0
      ZoW=0.0005
      Wstar=0.04*delw
      Tstar=0.04*delt
      Qstar=0.04*delq
      TVstar=Tstar*(1.0+0.61*qa)+0.61*ta_k*Qstar

!     Compute Richardson number.
      ri=g*zw*(delt+0.61*ta_k*delq)/(ta_k*delw*delw)

!     Fairall computes turbulent fluxes only if Ri< 0.25
      if ( ri .le. 0.25) then
!        Iterate until convergence or when IER is negative.  It usually
!        converges within four iterations.
         do iter=1,itermax
            if ( ier .ge. 0 ) then
!              Compute Monin-Obukhov stability parameter, Z/L.
               oL=g*kappa*TVstar/(ta_k*(1.0+0.61*qa)*Wstar*Wstar)
               ZWoL=zw*oL
               ZToL=zt*oL
               ZQoL=zq*oL

!              Evaluate stability functions at Z/L.
               wpsi=psi(1,ZWoL)
               tpsi=psi(2,ZToL)
               qpsi=psi(2,ZQoL)

!              Compute wind scaling parameters, Wstar.
               ZoW=0.011*Wstar*Wstar/g+0.11*vis_air/Wstar
               Wstar=delw*kappa/(log(zw/ZoW)-wpsi)

!              Computes roughness Reynolds number for wind (Rr), heat (Rt),
!              and moisture (Rq). Use Liu et al. (1976) look-up table to
!              compute "Rt" and "Rq" as function of "Rr".
               rr=ZoW*Wstar/vis_air
               if ((rr .ge. 0.0).and.(rr .lt. 1000.0)) then
                  do k=1,8
                     if ((liu_rr(k).le.rr).and.(rr .lt. liu_rr(k+1))) then
                        rt=liu_a(k,1)*rr**liu_b(k,1)
                        rq=liu_a(k,2)*rr**liu_b(k,2)
                     end if
                  end do

!                Compute heat and moisture scaling parameters,
!                Tstar and Qstar.
                  cff=vis_air/Wstar
                  ZoT=rt*cff
                  ZoQ=rq*cff
                  cff=kappa*fdg
                  Tstar=(delt)*cff/(log(zt/ZoT)-tpsi)
                  Qstar=(delq)*cff/(log(zq/ZoQ)-qpsi)

!                 Compute gustiness in wind speed.
                  TVstar=Tstar*(1.0+0.61*qa)+0.61*ta_k*Qstar
                  bf=-g/ta_k*Wstar*TVstar
                  if (bf .gt. 0) then
                     wgus=beta*(bf*Zabl)**r3
                  else
                     wgus=_ZERO_
                  end if
                  delw=sqrt(w*w+wgus*wgus)
               else
                  ier = -2
               end if
            end if
         end do

!        Compute transfer coefficients for momentun (Cd), heat (Ch),
!        and moisture (Ce).
         if (ier .ge. 0.0) then
            Wspeed=sqrt(w*w+wgus*wgus)
            Cd=Wstar*Wstar/(Wspeed*Wspeed)

!           Compute turbulent sensible heat flux (W/m2), qe.
!           out of ocean is negative
            qe=cpa*rhoa*Wstar*Tstar

!           compute sensible heatflux due to rain fall
            if (rain_impact) then
!              units of qs and qa - should be kg/kg
               rainfall=precip * 1000. ! (convert from m/s to kg/m2/s)
               x1 = 2.11e-5*(ta_k/kelvin)**1.94
               x2 = 0.02411*(1.0+ta*(3.309e-3-1.44e-6*ta))/(rhoa*cpa)
               x3 = qa * L /(rgas * ta_K * ta_K)
               cd_rain = 1.0/(1.0+const06*(x3*L*x1)/(cpa*x2))
               cd_rain = cd_rain*cpw*((tw-ta) + (qs-qa)*L/cpa)
               qe = qe - rainfall * cd_rain
            end if

!           Compute turbulent latent heat flux (W/m2), qh.
            qh=L*rhoa*Wstar*Qstar

!           Compute Webb correction (Webb effect) to latent heat flux
            upvel=-1.61*Wstar*Qstar-(1.0+1.61*qa)*Wstar*Tstar/ta_k
            qh=qh-rhoa*L*upvel*qa

!           calculation of evaporation/condensation in m/s
            if (rain_impact .and. calc_evaporation) then
               evap = rhoa/rho_0*Wstar*Qstar
            else
               evap = _ZERO_
            end if

!           Compute wind stress components (N/m2), Tau.
            cff=rhoa*Cd*Wspeed
            taux=(cff*u10)
            tauy=(cff*v10)

!           Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
            if ( rain_impact ) then
               tmp  = 0.85 * rainfall * Wspeed
               x=u10
               taux  = taux + tmp * sign(_ONE_,x)
               x=v10
               tauy  = tauy + tmp * sign(_ONE_,x)
            end if

         end if ! ier >0
      end if ! Ri < 0.25
   end if  !delw != 0.0

   return
   end subroutine fairall

   function psi(iflag, ZoL)
!=======================================================================
!                                                                      !
!  This function evaluates the stability function, PSI, for wind       !
!  speed (iflag=1) or for air temperature and moisture (iflag=2)       !
!  profiles as function of the stability parameter, ZoL (z/L).         !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!=======================================================================
!
      REALTYPE  :: psi
!
!  Imported variable declarations.
!
   integer,  intent(in) :: iflag
   REALTYPE, intent(in) :: ZoL
!
!  Local variable declarations.
!
   REALTYPE, parameter :: r3 = 1.0/3.0
   REALTYPE, parameter :: sqr3 = 1.7320508
   REALTYPE, parameter :: pi=3.141592653589
   REALTYPE            :: Fw, chic, chik, psic, psik

!  Initialize for the zero "ZoL" case.
!
   psi=0.0
!
!  Unstable conditions.
!
   if (ZoL .lt. 0.0) then
      chik=(1.0-16.0*ZoL)**0.25
      if (iflag .eq. 1) then
         psik=2.0*LOG(0.5*(1.0+chik))+LOG(0.5*(1.0+chik*chik))-   &
              2.0*ATAN(chik)+ 0.5*pi
      else if (iflag .eq. 2) then
            psik=2.0*LOG(0.5*(1.0+chik*chik))
      end if
!
!  For very unstable conditions, use free-convection (Fairall).
!
      chic=(1.0-12.87*ZoL)**r3
      psic=1.5*LOG(r3*(1.0+chic+chic*chic))-                    &
         sqr3*ATAN((1.0+2.0*chic)/sqr3)+ pi/sqr3
!
!  Match Kansas and free-convection forms with weighting Fw.
!
      Fw=1.0/(1.0+ZoL*ZoL)
      psi=Fw*psik+(1.0-Fw)*psic
!
!  Stable conditions.
!
   else if (ZoL .gt. 0.0) then
      psi=-4.7*ZoL
   end if

   return
   end function psi

!EOC
!-----------------------------------------------------------------------
!Copyright (C) 2007 - Adolf Stips
!-----------------------------------------------------------------------


!$Id$
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the short--wave radiation \label{sec:swr}
!
! !INTERFACE:
   subroutine short_wave_radiation(jul,secs,dlon,dlat,cloud,swr)
!
! !DESCRIPTION:
!  This subroutine calculates the short--wave net radiation based on
!  latitude, longitude, time, fractional cloud cover and albedo.
!  The albedo monthly values from \cite{Payne72} are given  here
!  as means of the values between
!  at 30$^{\circ}$ N and 40$^{\circ}$ N for the Atlantic Ocean
!  (hence the same latitudinal band of the Mediterranean Sea).
!  The basic formula for the short-wave radiation at the surface, $Q_s$,
!  has been taken from \cite{RosatiMiyacoda88}, who adapted the work
!  of \cite{Reed77} and \cite{SimpsonPaulson99}:
!
!  \begin{equation}
!  Q_s=Q_{tot} (1-0.62 C + 0.0019 \beta) (1-\alpha),
!  \end{equation}
!
!  with the total radiation reaching the surface under clear skies,
!  $Q_{tot}$, the fractional cloud cover, $C$, the solar noon altitude,
!  $\beta$, and the albedo, $\alpha$.
!  This piece of code has been taken the MOM-I (Modular Ocean Model)
!  version at the INGV (Istituto Nazionale di Geofisica e Vulcanologia,
!  see {\tt http://www.bo.ingv.it/}).
!
! !USES:
   use time, only: calendar_date
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
   REALTYPE, intent(in)                :: dlon,dlat
   REALTYPE, intent(in)                :: cloud
!
! !OUTPUT PARAMETERS:
   REALTYPE, optional, intent(out)     :: swr
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  $Log: short_wave_radiation.F90,v $
!  Revision 1.1  2007-09-25 10:06:10  kbk
!  modularized the airsea module - added Fairall method
!
!
! !LOCAL VARIABLES:
   REALTYPE, parameter       :: pi=3.14159265358979323846
   REALTYPE, parameter       :: deg2rad=pi/180.
   REALTYPE, parameter       :: rad2deg=180./pi

   REALTYPE                  :: solar=1350.
   REALTYPE                  :: eclips=23.439*deg2rad
   REALTYPE                  :: tau=0.7
   REALTYPE                  :: aozone=0.09

   REALTYPE                  :: th0,th02,th03,sundec
   REALTYPE                  :: thsun,coszen,zen,dzen,sunbet
   REALTYPE                  :: qatten,qzer,qdir,qdiff,qtot,qshort
   REALTYPE                  :: albedo
   integer                   :: jab
   integer                   :: yy,mm,dd
   REALTYPE                  :: yrdays,days,hour,tjul
   REALTYPE                  :: rlon,rlat

   integer                   :: yday(12) = &
                 (/ 0,31,59,90,120,151,181,212,243,273,304,334 /)

   REALTYPE                  :: alb1(20) = &
                 (/.719,.656,.603,.480,.385,.300,.250,.193,.164, &
                   .131,.103,.084,.071,.061,.054,.039,.036,.032,.031,.030 /)

   REALTYPE                  :: za(20) = &
                 (/90.,88.,86.,84.,82.,80.,78.,76.,74.,70.,  &
                   66.,62.,58.,54.,50.,40.,30.,20.,10.,0.0 /)

   REALTYPE                  :: dza(19)
   data           dza/8*2.0, 6*4.0, 5*10.0/
!EOP
!-----------------------------------------------------------------------
!BOC
!  from now on everything in radians
   rlon = deg2rad*dlon
   rlat = deg2rad*dlat

!  number of days in a year :
   call calendar_date(jul,yy,mm,dd)
   days=float(yday(mm))+float(dd)
   hour=1.0*secs/3600.
!kbk   if (mod(yy,4) .eq. 0 ! leap year I forgot
   yrdays=365.

   th0 = 2.*pi*days/yrdays
   th02 = 2.*th0
   th03 = 3.*th0
!  sun declination :
   sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0)         &
           - 0.006758*cos(th02) + 0.000907*sin(th02)                 &
           - 0.002697*cos(th03) + 0.001480*sin(th03)

!  sun hour angle :
   thsun = (hour-12.)*15.*deg2rad + rlon

!  cosine of the solar zenith angle :
   coszen =sin(rlat)*sin(sundec)+cos(rlat)*cos(sundec)*cos(thsun)
   if (coszen .le. 0.0) then
      coszen = 0.0
      qatten = 0.0
   else
      qatten = tau**(1./coszen)
   end if
   qzer  = coszen * solar
   qdir  = qzer * qatten
   qdiff = ((1.-aozone)*qzer - qdir) * 0.5
   qtot  =  qdir + qdiff

   tjul = (days-81.)/yrdays*2.*pi

!  sin of the solar noon altitude in radians :
   sunbet=sin(rlat)*sin(eclips*sin(tjul))+cos(rlat)*cos(eclips*sin(tjul))
!  solar noon altitude in degrees :
   sunbet = asin(sunbet)*rad2deg

!  calculates the albedo as a function of the solar zenith angle :
!  (after Payne jas 1972)
!  solar zenith angle in degrees :
   zen=(180./pi)*acos(coszen)
   if(zen .ge. 74.)then
      jab=.5*(90.-zen)+1.
   else if (zen .ge. 50.) then
      jab=.23*(74.-zen)+9.
   else
      jab=.10*(50.-zen)+15.
   endif

   dzen=(za(jab)-zen)/dza(jab)
   albedo=alb1(jab)+dzen*(alb1(jab+1)-alb1(jab))

!  radiation as from Reed(1977), Simpson and Paulson(1979)
!  calculates SHORT WAVE FLUX ( watt/m*m )
!  Rosati,Miyakoda 1988 ; eq. 3.8
!  clouds from COADS perpetual data set
#if 1
   qshort  = qtot*(1-0.62*cloud + .0019*sunbet)*(1.-albedo)
   if(qshort .gt. qtot ) then
      qshort  = qtot
   end if
#else
!  original implementation
   if(cloud .lt. 0.3) then
      qshort  = qtot
   else
      qshort  = qtot*(1-0.62*cloud + 0.0019*sunbet)*(1.-albedo)
   endif
#endif

   swr = qshort

   return
   end subroutine short_wave_radiation
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

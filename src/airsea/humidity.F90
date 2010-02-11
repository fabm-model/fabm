!$Id: humidity.F90,v 1.2 2007-10-02 10:14:08 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the humidity \label{sec:humidity}
!
! !INTERFACE:
   subroutine humidity(hum_method,hum,airp,tw,ta)
!
! !DESCRIPTION:
!
! This routine calculated the saturation vapour pressure at SST and at
! air temperature, as well as the saturation specific humidty and the 
! specific humidity. For the latter, four methods are implemented,
! and the method has to be chosen in the namelist file {\tt airsea.nml}
! as parameter {\tt hum\_method}, see \sect{sec:init-air-sea} for details. 
!
!
! !USES:
   use airsea_variables, only: kelvin,const06,rgas
   use airsea_variables, only: es,ea,qs,qa,rhoa
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: hum_method
   REALTYPE, intent(in)                :: hum,airp,tw,ta
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Adolf Stips, Hans Burchard & Karsten Bolding
!
! !DEFINED PARAMETERS:
   REALTYPE, parameter       :: a1=6.107799961
   REALTYPE, parameter       :: a2=4.436518521e-1
   REALTYPE, parameter       :: a3=1.428945805e-2
   REALTYPE, parameter       :: a4=2.650648471e-4
   REALTYPE, parameter       :: a5=3.031240396e-6
   REALTYPE, parameter       :: a6=2.034080948e-8
   REALTYPE, parameter       :: a7=6.136820929e-11
!
! !LOCAL VARIABLES:
   REALTYPE        :: rh,twet,twet_k,dew,dew_k
!EOP
!-----------------------------------------------------------------------
!BOC
!  saturation vapor pressure - using SST
   es = a1 +tw*(a2+tw*(a3+tw*(a4+tw*(a5+tw*(a6+tw*a7)))))
   es = es * 100.0 ! Conversion millibar --> Pascal

!  correction for seawater, following Kraus 1972
!  correcting for salt water assuming 98% RH
   es=0.98 * es
!  saturation specific humidity
   qs = const06*es/(airp-0.377*es)

!  must be also calcuated for airtemperature, depending on humidity input
!  see ../ncdf/ncdf_meteo.F90 for defined constants
   select case (hum_method)
      case (1) ! relative humidity in % given
         rh = 0.01 * hum
!        saturation vapor pressure at that air temperature
         ea = a1 +ta*(a2+ta*(a3+ta*(a4+ta*(a5+ta*(a6+ta*a7)))))
         ea = ea * 100.0 ! Conversion millibar --> Pascal
!        get actual vapor pressure
         ea = rh * ea
!        convert to specific humidity
         qa = const06*ea/(airp-0.377*ea)
      case (2)  ! Specific humidity from wet bulb temperature
!        calculate the SVP at wet bulb temp then
!        use the psychrometer formula to get vapour pressure
!        See Smithsonian Met tables 6th Edition pg 366 eqn 3
!        Make sure this is in degC
         if (hum .lt. 100 ) then
            twet_k=hum + kelvin
            twet=hum
         else
            twet=hum - kelvin
            twet_k=hum
         end if
!        saturation vapor pressure at wet bulb temperature
         ea = a1 +twet*(a2+twet*(a3+twet*(a4+twet*(a5+twet*(a6+twet*a7)))))
         ea = ea * 100.0 ! Conversion millibar --> Pascal
!        actual vapor pressure
         ea = ea - 6.6e-4*(1+1.15e-3*twet)*airp*(ta-twet)
!        specific humidity in kg/kg
         qa = const06*ea/(airp-0.377*ea)
      case (3)  ! Specific humidity from dew point temperature
         if (hum .lt. 100.) then
            dew = hum
            dew_k = hum + kelvin
         else
            dew = hum - kelvin
            dew_k = hum
         end if
         ea = a1 +dew*(a2+dew*(a3+dew*(a4+dew*(a5+dew*(a6+dew*a7)))))
         ea = ea * 100.0 ! Conversion millibar --> Pascal
         qa = const06*ea/(airp-0.377*ea)
      case (4)
!        specific humidity in kg/kg is given
         qa = hum
!        actual water vapor pressure in Pascal
         ea = qa *airp/(const06+0.378*qa)
      case default
         FATAL 'not a valid hum_method'
         stop 'bulk_fluxes()'
   end select

   rhoa = airp/(rgas*(ta+kelvin)*(1.0+const06*qa))

   return
   end subroutine humidity
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

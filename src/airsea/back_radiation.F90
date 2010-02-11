!$Id: back_radiation.F90,v 1.1 2007-09-25 10:06:10 kbk Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the long-wave back radiation \label{sec:back-rad}
!
! !INTERFACE:
   subroutine back_radiation(method,dlat,tw,ta,cloud,qb)
!
! !DESCRIPTION:
!
! Here, the long-wave back radiation is calculated by means of one out
! of four methods, which depend on the value given to the parameter
! {\tt method}: 
! {\tt method}=1: \cite{Clarketal74},  
! {\tt method}=2: \cite{HastenrathLamb78},  
! {\tt method}=3: \cite{Bignamietal95},  
! {\tt method}=4: \cite{BerliandBerliand52}.
! It should ne noted that the latitude must here be given in degrees.
!
! !USES:
   use airsea_variables, only: emiss,bolz
   use airsea_variables, only: ea,qa
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: method
   REALTYPE, intent(in)                :: dlat,tw,ta,cloud
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: qb
!
! !REVISION HISTORY:
!  Original author(s): Adols Stips, Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:

   integer, parameter   :: clark=1      ! Clark et al, 1974
   integer, parameter   :: hastenrath=2 ! Hastenrath and Lamb, 1978
   integer, parameter   :: bignami=3    ! Bignami et al., 1995 - Medsea
   integer, parameter   :: berliand=4   ! Berliand and Berliand, 1952 - ROMS


   REALTYPE, parameter, dimension(91)  :: cloud_correction_factor = (/ &
     0.497202,     0.501885,     0.506568,     0.511250,     0.515933, &
     0.520616,     0.525299,     0.529982,     0.534665,     0.539348, &
     0.544031,     0.548714,     0.553397,     0.558080,     0.562763, &
     0.567446,     0.572129,     0.576812,     0.581495,     0.586178, &
     0.590861,     0.595544,     0.600227,     0.604910,     0.609593, &
     0.614276,     0.618959,     0.623641,     0.628324,     0.633007, &
     0.637690,     0.642373,     0.647056,     0.651739,     0.656422, &
     0.661105,     0.665788,     0.670471,     0.675154,     0.679837, &
     0.684520,     0.689203,     0.693886,     0.698569,     0.703252, &
     0.707935,     0.712618,     0.717301,     0.721984,     0.726667, &
     0.731350,     0.736032,     0.740715,     0.745398,     0.750081, &
     0.754764,     0.759447,     0.764130,     0.768813,     0.773496, &
     0.778179,     0.782862,     0.787545,     0.792228,     0.796911, &
     0.801594,     0.806277,     0.810960,     0.815643,     0.820326, &
     0.825009,     0.829692,     0.834375,     0.839058,     0.843741, &
     0.848423,     0.853106,     0.857789,     0.862472,     0.867155, &
     0.871838,     0.876521,     0.881204,     0.885887,     0.890570, &
     0.895253,     0.899936,     0.904619,     0.909302,     0.913985, &
     0.918668 /)

   REALTYPE                  :: ccf
   REALTYPE                  :: x1,x2,x3
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  calculate cloud correction factor,fortran counts from 1 !
   ccf= cloud_correction_factor(nint(abs(dlat))+1)

   select case(method)
      case(clark)
!        Clark et al. (1974) formula.
!        unit of ea is Pascal, must hPa
!        Black body defect term, clouds, water vapor correction
         x1=(1.0-ccf*cloud*cloud)*(tw**4)
         x2=(0.39-0.05*sqrt(ea*0.01))
!        temperature jump term
         x3=4.0*(tw**3)*(tw-ta)
         qb=-emiss*bolz*(x1*x2+x3)
      case(hastenrath) ! qa in g(water)/kg(wet air)
!        Hastenrath and Lamb (1978) formula.
         x1=(1.0-ccf*cloud*cloud)*(tw**4)
         x2=(0.39-0.056*sqrt(1000.0*qa))
         x3=4.0*(tw**3)*(tw-ta)
         qb=-emiss*bolz*(x1*x2+x3)
      case(bignami)
!        Bignami et al. (1995) formula (Med Sea).
!        unit of ea is Pascal, must hPa
         ccf=0.1762
         x1=(1.0+ccf*cloud*cloud)*ta**4
         x2=(0.653+0.00535*(ea*0.01))
         x3= emiss*(tw**4)
         qb=-bolz*(-x1*x2+x3)
      case(berliand)
!        Berliand & Berliand (1952) formula (ROMS).
         x1=(1.0-0.6823*cloud*cloud)*ta**4
         x2=(0.39-0.05*sqrt(0.01*ea))
         x3=4.0*ta**3*(tw-ta)
         qb=-emiss*bolz*(x1*x2+x3)
      case default
   end select

   return
   end subroutine back_radiation
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

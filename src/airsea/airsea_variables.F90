!$Id: airsea_variables.F90,v 1.4 2008-05-02 11:49:12 kb Exp $
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: airsea_variables\label{sec:airsea-variables}
!
! !INTERFACE:
   module airsea_variables
!
! !DESCRIPTION:
!
!  Here, number of public variables in the airsea module is declared. 
!
! !USES:
   IMPLICIT NONE

!  default: all is private.
   private
!
! !PUBLIC DATA MEMBERS:
   REALTYPE, public, parameter         :: cpa=1008.
   REALTYPE, public, parameter         :: cpw=3985.
   REALTYPE, public, parameter         :: emiss=0.97
   REALTYPE, public, parameter         :: bolz=5.67e-8
   REALTYPE, public, parameter         :: kelvin=273.16
   REALTYPE, public, parameter         :: const06=0.62198
   REALTYPE, public, parameter         :: rgas = 287.1    ! 
   REALTYPE, public, parameter         :: g = 9.81        ! [m/s2]
   REALTYPE, public, parameter         :: rho_0 = 1025.   ! [kg/m3]
   REALTYPE, public, parameter         :: kappa = 0.41    ! von Karman
   REALTYPE, public                    :: es=_ZERO_
   REALTYPE, public                    :: ea=_ZERO_
   REALTYPE, public                    :: qs=_ZERO_
   REALTYPE, public                    :: qa=_ZERO_
   REALTYPE, public                    :: L=_ZERO_
   REALTYPE, public                    :: rhoa=_ZERO_
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding, Hans Burchard
!
!  $Log: airsea_variables.F90,v $
!  Revision 1.4  2008-05-02 11:49:12  kb
!  explicitely initialise variable to 0 - solves Mac bug
!
!  Revision 1.3  2007-12-21 12:38:03  kb
!  added precip/evap to kondo + cleaned
!
!  Revision 1.2  2007-10-02 10:14:08  kbk
!  fixed rhoa calculation - rgas in airsea_variables module
!
!  Revision 1.1  2007-09-25 10:06:10  kbk
!  modularized the airsea module - added Fairall method
!
!
!EOP
!-----------------------------------------------------------------------

   end module airsea_variables

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

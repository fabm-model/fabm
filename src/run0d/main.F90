#include<cppdefs.h>
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: 0D shell based on GOTM --- the main program  \label{sec:main}
!
! !INTERFACE:
   program main
!
! !DESCRIPTION:
! TODO
!
! !USES:
   use time
   use run0d
!
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   character(LEN=8)          :: datestr
   real                      :: t1=-1,t2=-1
!
!-----------------------------------------------------------------------
!BOC
#ifdef FORTRAN95
   call CPU_Time(t1)
#endif
   call Date_And_Time(datestr,timestr)
   STDERR LINE
   STDERR '0D shell based on GOTM ver. ',RELEASE,': Started on  ',datestr,' ',timestr
   STDERR LINE

   call init_run()
   call time_loop()
   call clean_up()

#ifdef FORTRAN95
   call CPU_Time(t2)
#endif
   call Date_And_Time(datestr,timestr)
   STDERR LINE
   STDERR '0D shell based on GOTM ver. ',RELEASE,': Finished on ',datestr,' ',timestr
#ifdef FORTRAN95
   STDERR 'CPU-time was in loop:  ',t2-t1,' seconds'
   STDERR 'Sim-time/CPU-time:     ',simtime/(t2-t1)
#endif
   STDERR LINE
   STDERR LINE

   end
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

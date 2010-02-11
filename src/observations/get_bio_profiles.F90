!$Id: get_bio_profiles.F90,v 1.3 2009-03-19 09:36:32 kb Exp $
#ifdef BIO
#include "cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_bio_profiles
!
! !INTERFACE:
   subroutine get_bio_profiles(unit,jul,secs,nlev,z)
!
! !DESCRIPTION:
!  This routine is responsible for providing sane values to `observed'
!  biological profiles.
!  The subroutine is called in the {\tt get\_all\_obs{}} subroutine
!  as part of the main integration loop.
!  In case of observations from file the temporal interpolation is
!  done in this routine.
!
! !USES:
   use time
   use observations, only: init_saved_vars,read_profiles,bioprofs
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in):: unit
   integer, intent(in):: jul,secs
   integer, intent(in):: nlev
   REALTYPE, intent(in):: z(0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!  $Log: get_bio_profiles.F90,v $
!  Revision 1.3  2009-03-19 09:36:32  kb
!  also work with BIO=false
!
!  Revision 1.2  2007-09-11 13:24:32  jorn
!  added stop after fatal error reading profile
!
!  Revision 1.1  2007-06-19 10:38:03  kbk
!  initialise biological profiles from external file
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: rc,n
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,dt
   integer, save             :: jul1,secs1
   integer, save             :: jul2=0,secs2=0
   integer, save             :: cols
   integer, save             :: lines=0
   integer, save             :: nprofiles=0
   logical, save             :: one_profile=.false.
   REALTYPE, save, dimension(:,:), allocatable :: prof1,prof2,alpha
!
!-----------------------------------------------------------------------
!BOC
   if (init_saved_vars) then
      jul2=0
      secs2=0
      cols=size(bioprofs,1)
      lines=0
      nprofiles=0
      one_profile=.false.
   end if

   if ( .not. allocated(prof1)) then
      allocate(prof1(0:nlev,cols),stat=rc)
      if (rc /= 0) stop 'get_bio_profiles: Error allocating memory (prof1)'
      prof1 = 0.
   end if
   if ( .not. allocated(prof2)) then
      allocate(prof2(0:nlev,cols),stat=rc)
      if (rc /= 0) stop 'get_bio_profiles: Error allocating memory (prof2)'
      prof2 = 0.
   end if
   if ( .not. allocated(alpha)) then
      allocate(alpha(0:nlev,cols),stat=rc)
      if (rc /= 0) stop 'get_bio_profiles: Error allocating memory (alpha)'
   end if

!  This part initialises and reads in new values if necessary.
   if(.not. one_profile .and. time_diff(jul2,secs2,jul,secs) .lt. 0) then
      do
         jul1 = jul2
         secs1 = secs2
         prof1 = prof2
         call read_profiles(unit,nlev,cols,yy,mm,dd,hh,min,ss,z,prof2,lines,rc)
         if(rc .ne. 0) then
            if(nprofiles .eq. 1) then
               LEVEL3 'Only one set of biological profiles are present.'
               one_profile = .true.
               do n=1,cols
                  bioprofs(n,:) = prof1(:,n)
               end do
            else
               FATAL 'Error reading biological profiles around line #',lines
               stop 'get_bio_profiles'
            end if
            EXIT
         else
            nprofiles = nprofiles + 1
            call julian_day(yy,mm,dd,jul2)
            secs2 = hh*3600 + min*60 + ss
            if(time_diff(jul2,secs2,jul,secs) .gt. 0) EXIT
         end if
      end do
      if( .not. one_profile) then
         dt = time_diff(jul2,secs2,jul1,secs1)
         alpha = (prof2-prof1)/dt
      end if
   end if

!  Do the time interpolation - only if more than one profile
   if( .not. one_profile) then
      t  = time_diff(jul,secs,jul1,secs1)
      do n=1,cols
         bioprofs(n,:) = prof1(:,n) + t*alpha(:,n)
      end do
   end if

   return
   end subroutine get_bio_profiles
!EOC
#endif

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

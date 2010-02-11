!$Id: bio_cl.F90,v 1.2 2009-12-09 13:11:54 kb Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_cl --- Chlorination model \label{sec:bio-cl}
!
! !INTERFACE:
   module bio_cl
!
! !DESCRIPTION:
!  Chlorination
!
! !USES:
!  default: all is private.
   use bio_var
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_cl, init_var_cl, do_bio_cl, end_bio_cl
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Hannes Rennau and Karsten Bolding
!                      with input from John Aldrige
!
!  $Log: bio_cl.F90,v $
!  Revision 1.2  2009-12-09 13:11:54  kb
!  added GUI for BIO/CL model - Rennau
!
!  Revision 1.1  2009-11-11 13:08:55  kb
!  added chlorination model - Rennau
!
!
! !LOCAL VARIABLES:
!  from a namelist
   integer                   :: cl_method
   REALTYPE                  :: cl_initial=0.1
   REALTYPE                  :: s_initial=0.58
   REALTYPE                  :: k1=0.000042
   REALTYPE                  :: k2=0.013
   REALTYPE                  :: wk1=-0.5757
   REALTYPE                  :: wk2=1.

   REALTYPE, public          :: cl0=0.
   REALTYPE                  :: s0=0.
   REALTYPE                  :: w_cl=-0.
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the cl module
!
! !INTERFACE:
   subroutine init_bio_cl(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio\_cl.nml} is read and
!  various variables (rates and settling velocities)
!  are transformed into SI units.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: namlst
   character(len=*), intent(in)   :: fname
   integer,          intent(in)   :: unit
!
! !REVISION HISTORY:
!  Original author(s): Hannes Rennau

! !LOCAL VARIABLES:
   namelist /bio_cl_nml/ numc, &
                      numc,cl_method,cl_initial,s_initial,k1,k2, &
                      wk1,wk2,cl0,s0,w_cl
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_cl'

   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_cl_nml,err=99)
   close(namlst)
   LEVEL2 'read_test'
   LEVEL3 'namelist "', fname, '" read'
 
   select case (cl_method)
      case (1) ! Cefas
         numc=2
      case (2) ! Wang
         numc=1
      case default
         stop 'init_bio_cl(): non valid cl_method'
   end select
 
!  Conversion from day to second
   w_cl = w_cl  /secs_pr_day

!  initialize variable descriptions

   call bio_alloc_info

   var_names(1) = 'cl'
   var_units(1) = 'mg/l'
   var_long(1) = 'chlorine concentration'

   if ( cl_method .eq. 1 ) then
      var_names(2) = 'cl_rc'
      var_units(2) = 'mg/l'
      var_long(2) = 'reactive compound (chlorine)'
   end if

   LEVEL3 'module initialized'

   return

98 LEVEL2 'I could not open bio_cl.nml'
   LEVEL2 'If thats not what you want you have to supply bio_cl.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_cl.nml'
   return
99 FATAL 'I could not read bio_cl.nml'
   stop 'init_bio_cl'
   end subroutine init_bio_cl
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_cl
!
! !DESCRIPTION:
!  Here, the the initial conditions are set and the settling velocities are
!  transferred to all vertical levels. All concentrations are declared
!  as non-negative variables, and it is defined which variables would be
!  taken up by benthic filter feeders.
!
! !USES:
   IMPLICIT NONE

! !REVISION HISTORY:
!  Original author(s):  Hannes Rennau

! !LOCAL VARIABLES:
  integer                    :: i
!EOP
!-----------------------------------------------------------------------
!BOC

   cc(1,:) = cl_initial
   ws(1,:) = w_cl
   posconc(1) = 1

   if ( cl_method .eq. 1 ) then
      cc(2,:) = s_initial
      ws(2,:) = w_cl
      posconc(2) = 1
   end if

   LEVEL3 'CL variables initialised ...'

   return

   end subroutine init_var_cl
!EOC

#if 0
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: decay of chlorine is calculated
!            waiting for further information on process...
! !INTERFACE:
   REALTYPE function dcl(cl,s,k1,k2)
!
! !DESCRIPTION:
! Here, the ...
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: cl,s,k1,k2
!
! !REVISION HISTORY:
!  Original author(s):  Hannes Rennau
!
!EOP
!-----------------------------------------------------------------------
!BOC

   dcl=k1*cl**2*s+k2*cl
   return
   end function dcl
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: decay of chlorine is calculated (Wang et al. 2008)
!
! !INTERFACE:
   REALTYPE function dcl_wang(cl,wk1,wk2)
!
! !DESCRIPTION:
! Here, the ...
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: cl,wk1,wk2
!
! !REVISION HISTORY:
!  Original author(s):  Hannes Rennau
!
!EOP
!-----------------------------------------------------------------------
!BOC

   dcl_wang = cl-(cl*2.7183**(wk1/3600.))/wk2
   return
   end function dcl_wang
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: decay of chlorine reactive compound is calculated
!            waiting for further information on process...
! !INTERFACE:
   REALTYPE function ds(cl,s,k1)
!
! !DESCRIPTION:
! Here, the ...
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: cl,s,k1
!
! !REVISION HISTORY:
!  Original author(s):  Hannes Rennau
!
!EOP
!-----------------------------------------------------------------------
!BOC

   ds=k1*cl**2*s
   return
   end function ds
!EOC

#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of CL model
!
! !INTERFACE:
   subroutine do_bio_cl(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: numc,nlev
   REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hannes Rennau
!
! !LOCAL VARIABLES:
   integer                    :: ci
   REALTYPE                   :: x
!EOP
!-----------------------------------------------------------------------
!BOC
   pp = _ZERO_

   select case (cl_method)
      case (1) ! Cefas
         do ci=1,nlev
!           dcl=k2*cl**2*s+k1*cl
            x=k2*cc(1,ci)*cc(1,ci)*cc(2,ci)
            dd(1,1,ci)=x+k1*cc(1,ci)
            dd(1,2,ci) = _ZERO_
!           ds=k2*cl**2*s
            dd(2,2,ci)=x
         end do
      case (2) ! Wang
         do ci=1,nlev
!           dcl_wang = cl-(cl*2.7183**(wk1/3600.))/wk2
            dd(1,1,ci)=cc(1,ci)-(cc(1,ci)*2.7183**(wk1/3600.))/wk2
         end do
    end select

   return
   end subroutine do_bio_cl
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_cl
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Hannes Rennau
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_bio_cl
!EOC

!-----------------------------------------------------------------------

   end module bio_cl

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

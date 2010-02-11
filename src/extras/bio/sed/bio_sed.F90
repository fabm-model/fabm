!$Id: bio_sed.F90,v 1.2 2008-07-08 10:04:23 lars Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_sed --- simple suspended matter model \label{sec:biosed}
!
! !INTERFACE:
   module bio_sed
!
! !DESCRIPTION:
!  This is a simple suspended matter model with one non-dimensional
!  state variable called {\tt conc}. The suspended matter is subject
!  to a constant settling velocity, has no surface fluxes of suspended matter,
!  but the suspended matter may be taken out at the bed, if the mussel
!  module of GOTM is activated. No right-hand side process terms are 
!  involved here. Note that this module has an Eulerian version, 
!  corresponding to a simple advection-diffusion equation with constant
!  settling velocity, and a Langrangian particle version, in which
!  particle diffusion is performed by random-walking the particles 
!  as described in \sect{sec:lagrangian}.
!
! !USES:
!  default: all is private.
   use bio_var
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_sed, init_var_sed, end_bio_sed
   public do_bio_sed_eul, do_bio_sed_par
!
! !REVISION HISTORY:!
!  Original author(s): Hans Burchard, Lars Umlauf, Karsten Bolding
!
!  $Log: bio_sed.F90,v $
!  Revision 1.2  2008-07-08 10:04:23  lars
!  changed initialization and particle support
!
!  Revision 1.1  2008-03-26 08:51:44  kb
!  new directory based bio structure
!
!  Revision 1.9  2007-04-18 07:35:26  kbk
!  to avoid F95 warning
!
!  Revision 1.8  2007-01-06 11:49:15  kbk
!  namelist file extension changed .inp --> .nml
!
!  Revision 1.7  2006-10-26 13:12:46  kbk
!  updated bio models to new ode_solver
!
!  Revision 1.6  2005-12-02 20:57:27  hb
!  Documentation updated and some bugs fixed
!
!  Revision 1.5  2005-11-17 09:58:18  hb
!  explicit argument for positive definite variables in diff_center()
!
!  Revision 1.4  2005/09/19 21:03:31  hb
!  pp and dd properly set to zero
!
!  Revision 1.3  2004/08/02 08:34:36  hb
!  updated init routines to reflect new internal bio interface
!
!  Revision 1.2  2004/07/30 09:22:20  hb
!  use bio_var in specific bio models - simpliefied internal interface
!
!  Revision 1.1  2003/10/28 10:22:45  hb
!  added support for sedimentation only 1 compartment bio model
!
!EOP
!-----------------------------------------------------------------------!
!
! !LOCAL VARIABLES:
!  from a namelist
   REALTYPE                  :: C_initial=4.5
   REALTYPE                  :: w_C=-5.787037e-05

! field IDs 
  integer, parameter                          ::  LoadInd=1


   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio_sed(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio\_sed.nml} (mainly including
!  settling velocity and initial value) is read
!  and the settling velocity is converted to SI units.
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
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------

! !LOCAL VARIABLES:
   namelist /bio_sed_nml/ numc,ntype,nprop,c_initial,w_C

!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_sed'

   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_sed_nml,err=99)
   close(namlst)

   LEVEL3 'namelist "', fname, '" read'

!  conversion from day to second
   w_C = w_C/secs_pr_day

!  initialize variable descriptions

   call bio_alloc_info

!  initialize field descriptions
   var_names(1) = 'conc'
   var_units(1) = '-'
   var_long(1)  = 'sediment concentration'


   LEVEL3 'module initialized'

   return

98 LEVEL2 'I could not open bio_sed.nml'
   LEVEL2 'If thats not what you want you have to supply bio_sed.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_sed.nml'
   return
99 FATAL 'I could not read bio_sed.nml'
   stop 'init_bio_sed'
   end subroutine init_bio_sed
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_sed
!
! !DESCRIPTION:
!  Here, fields names for the output variables are assigned
!  (in this case only one variable describing the sediment
!  concentration {\tt conc} exisits), and the concentration
!  field is initialized with a uniform concentration as 
!  specified in the namelist file {\tt bio\_sed.nml}. If the 
!  particle solver is used, all particles are distributed 
!  homogeneously over the water column.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf,  Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------!
! !LOCAL VARIABLES:

   integer                    :: nt

!-----------------------------------------------------------------------!
!BOC

! fill-in concentration arays and sinking speed
  if (bio_eulerian) then

     cc(1,:)    = C_initial
     
     ws(1,:)    = w_C
     
     posconc(1) = 1

  else
     
     !  all fields empty
     cc = _ZERO_
     
     !  constant sinking
     ws(1,:)    = w_C

     ! initalize particle positions
     call init_par_sed

     ! constant particle load
     nt = 1

     par_prop(:,LoadInd,nt) = (ztop-zbot)*C_initial/npar
     
  endif
  

#if 0
   mussels_inhale(1) = .true.
#endif

   LEVEL3 'variables initalized'

   return

   end subroutine init_var_sed
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fill initial particle distribution
!
! !INTERFACE:
   subroutine init_par_sed()
!
! !DESCRIPTION:
! Particles are distributed homogeneously over the whole water column. 
! Indices of vertical grid cells are assigend to all particles, and the
! particles are marked active.  
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf, Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   integer                   :: i,j,n

!-----------------------------------------------------------------------
!BOC


!  create homogeneous particle distribution
   do j=1,ntype
      do n=1,npar
         par_z(n,j) = zbot + n/float(npar+1)*(ztop-zbot)
      end do
   end do

!  assign cell indices to particles
   do j=1,ntype
      do n=1,npar
         do i=1,nlev
            if (zlev(i) .gt. par_z(n,j)) exit
         end do
         par_ind(n,j) = i
         par_act(n,j)  = .true.
      end do
   end do


   return
   end subroutine init_par_sed
!EOC
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of geobiochemical model
!
! !INTERFACE:
   subroutine do_bio_sed_eul(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! This routine sets the sinks and sources of this simple suspended
! matter module to zero, since sediment is considered as a passive
! tracer here.
!
! !USES:
   IMPLICIT NONE

! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: numc,nlev
   REALTYPE, intent(in)                 :: cc(1:numc,0:nlev)

! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  no right hand sides necessary

   pp=_ZERO_
   dd=_ZERO_

   return
   end subroutine do_bio_sed_eul
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update particle model
!
! !INTERFACE:
   subroutine do_bio_sed_par
!
! !DESCRIPTION:
!  Only provided for completeness since particle properties are 
! assumed passive here and do not change. 
!
! !USES:
   IMPLICIT NONE
!
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   call do_statistics
   
   return
   end subroutine do_bio_sed_par
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute bin statistics
!
! !INTERFACE:
   subroutine do_statistics
!
! !DESCRIPTION:
!  Computes some statistical properties of particles
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf & Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   integer                    :: ind,np,nt,i
!
!-----------------------------------------------------------------------
!BOC

   cc(1,:) = _ZERO_

   nt = 1;   ! only one particle type here

   do np=1,npar
      if (par_act(np,nt)) then

         ind = par_ind(np,nt) 

         ! add load concentrations
         cc(1,ind) = cc(1,ind) + par_prop(np,LoadInd,nt)

      end if

   end do

   ! compute volume averages
   do i=1,nlev
      if (cc(1,i) /= 0) then
         cc(1,i) = cc(1,i)/h(i)
      else
         cc(:,i) = _ZERO_
      end if
   end do
      
   return
   end subroutine do_statistics
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_sed
!
! !DESCRIPTION:
!  Nothing done at the moment
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_bio_sed

!EOC
!-----------------------------------------------------------------------

   end module bio_sed

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

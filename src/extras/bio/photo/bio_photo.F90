!$Id$
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_photo --- simple photoinhibition model \label{sec:bio-photo}
!
! !INTERFACE:
   module bio_photo
!
! !DESCRIPTION:
! This module calculates the inhibition parameter,
! $Y$, and the photosynthetic production rates $P$ according to a
! theory of \cite{DenmanMarra1986a} as adapted by
! \cite{Nagaietal2003a}.  The photosynthetically available radiation
! $PAR$ is given as \begin{equation} PAR(z)=I_0 (1-a)
! \exp\left(\frac{z}{\gamma_2} \right) \; , \end{equation} where $z$
! denotes the vertical coordinate, $(1-a)$ the short-wave fraction of
! visible light identified with $PAR$, $\gamma_2$ the penetration
! depth of the short-wave fraction, and $I_0$ the incoming solar
! radiation. The parameters $a$ and $\gamma_2$ have to be provided
! separately here via the namelist of this model, but they should be
! made consistent with the values chosen in {\tt obs.inp} via the 
! Jerlov water class.
!
! Note that radiation units used by \cite{Nagaietal2003a}
! have been converted to units W m$^{-2}$ s$^{-1}$ using the
! conversion relation $1$ W m$^{-2}$ s$^{-1}$ $ \approx 0.2174 \, \mu$E
! m$^{-2}$ s$^{-1}$, which is strictly valid only for a wave length of 
! $550$ nm. 
!
! Instantaneous photosynthetic production rates for fully uninhibited and
! fully inhibited cells are
! \begin{equation}
! \begin{array}{l}
! \displaystyle
! P_d = P_{dm} \left(1-\exp\left(-\frac{PAR}{E_d}  \right)  \right) \; , \\ \\
! \displaystyle
! P_l = P_{lm} \left(1-\exp\left(-\frac{PAR}{E_l}  \right)  \right) \; ,
! \end{array}
! \end{equation}
! respectively, with the maximum production rates, $P_{dm}$ and
! $P_{lm}$, and the saturation values, $E_d$ and $E_l$. In contrast to 
! \cite{Nagaietal2003a}, here production is in carbon units (converted
! from the oxygen units by using the appropriate molar weights).
!
! The local production rate is then calculated by
! \begin{equation}
!   P=P_d + Y (P_l-P_d) \; ,
! \end{equation}
! where $0 \le Y \le 1$ denotes the instantaneous inhibition
! factor determining whether cells suffer from photoinhibition.
! The inhibition parameter $Y$ is computed from a linear first-order 
! differential equation,
! \begin{equation}
! \frac{\mbox{d} Y}{\mbox{d} t} = \frac{1}{\gamma} (X-Y) \; ,
! \end{equation}
! where $\gamma$ is the response time scale and $X$ is the 
! inhibition factor after full adaption to the local conditions. 
! It is computed from 
! \begin{equation}
!   X=1-\exp \left(-\left(\frac{\max\{PAR,E_b\}-E_b}{E_b}\right)^2 \right).
! \end{equation}
!
! From these production parameters, the cell carbon content is then computed
! according to 
! \begin{equation}
!  \deriv{C}{t} = P - \mu C \; ,
! \end{equation}
! where $\mu$ denotes a constant respiration rate.
!
! If the parameter {\tt photo\_mig} is set {\tt .true.} particles
! migrate vertically according to the following algorithm. From
! 06:00 am to 06:00 pm particles move upwards unless their individual
! photoinhibition parameter, $Y$, exceeds the threshold $Y_h$. Those
! particles move downward, i.e.\ away from bright regions. From 
! 06:00 pm to 06:00 am, all particles move downward. For all types
! of motion the migration speed, $w_c$, is identical and constant.
!
! In the present version, this model is implemented for individual
! Lagrangian particles ('cells'). Averages of the parameters $X$ and
! $Y$ are computed for each grid volume for output. Volume
! concentrations of $P$ and $C$ are computed for output by summing
! over all cell concentrations inside a grid volume, and dividing by
! the grid volume. Since the number of particles is arbitrary in these
! simulations, the concentration values are only determined within a
! constant factor.
!
! !USES:
!  default: all is private.
   use bio_var
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_photo, init_var_photo, end_bio_photo
   public do_bio_photo_eul,do_bio_photo_par
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf, Hans Burchard, Karsten Bolding
!
!  $Log: bio_photo.F90,v $
!  Revision 1.3  2008-10-17 08:48:29  lars
!  fixed bug in carbon time stepping
!
!  Revision 1.2  2008-07-08 10:40:28  lars
!  added external info string allocation
!
!  Revision 1.1  2008-07-08 10:17:25  lars
!  first version of particle model for photo response
!

! !LOCAL VARIABLES:
!
!  ______namelist variables______

!  do photomigration
   logical                                    ::  photo_mig=.false.

!  photoresponse time (s)
   REALTYPE                                   ::  gamma=     3600. 
  
!  initial carbon concentration (pg C/cell)
   REALTYPE                                   ::  cinit=     0.

!  uninhibited production (pg C/cell/hour)
   REALTYPE                                   ::  pdm=       4.0 

!  inhibited production (pg C/cell/hour)
   REALTYPE                                   ::  plm=       0.25

!  respiration rate (1/hour)
   REALTYPE                                   ::  mu=        0.1

!  Radiation treshhold (W/m2)
   REALTYPE                                   ::  ed=        163.05

!  Radiation treshhold (W/m2)
   REALTYPE                                   ::  el=        163.05

!  Inhibition treshhold (W/m2)
   REALTYPE                                   ::  eb=        43.48

!  Inhib. threshold for upward migration
   REALTYPE                                   ::  yl=        0.2

!  Inhib. threshold for downward migration
   REALTYPE                                   ::  yh=        0.2

!  Migration speed (m/s)
   REALTYPE                                   ::  w_c=       0.0003

!  (1-aa) is PAR fraction of visible light
   REALTYPE                                   ::  aa=        0.78   

!  penetration depth of PAR (m)
   REALTYPE                                   ::  g2=        7.9

!EOP
!-----------------------------------------------------------------------


!  field IDs 
   integer, parameter                          ::  XInd=1,YInd=2
   integer, parameter                          ::  ProdInd=3,CarbInd=4


   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise photoadaption
!
! !INTERFACE:
   subroutine init_bio_photo(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here, the namelist {\tt bio\_photo.nml} is read, some 
!  variables are converted to SI units, and the variable
!  descriptions (names, units) are set. 
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
!  Original author(s): Hans Burchard, Lars Umlauf, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:

   namelist /bio_photo_nml/ numc,ntype,nprop,photo_mig,gamma,           &
                            pdm,plm,cinit,mu,ed,el,eb,yl,yh,            &
                            w_c,aa,g2
!-----------------------------------------------------------------------                      
!BOC
   LEVEL2 'init_bio_photo'

   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_photo_nml,err=99)
   close(namlst)

   LEVEL3 'namelist "', fname, '" read'

!  conversion from day to second
   pdm  = pdm  /secs_pr_hour
   plm  = plm  /secs_pr_hour
   mu   = mu   /secs_pr_hour

!  initialize variable descriptions

   call bio_alloc_info


   var_names(1) = 'Np'
   var_units(1) = 'counts/volume'
   var_long(1)  = 'number of particles per volume'

   var_names(2) = 'X'
   var_units(2) = '-'
   var_long(2)  = 'fully adapted inhibition factor'

   var_names(3) = 'Y'
   var_units(3) = '-'
   var_long(3)  = 'inhibition factor'

   var_names(4) = 'Prod'
   var_units(4) = 'pg C/hour/vol'
   var_long(4)  = 'production'

   var_names(5) = 'Carb'
   var_units(5) = 'pg C/vol'
   var_long(5)  = 'carbon concentration'


   LEVEL3 'module initialized'

   return

98 LEVEL2 'I could not open "bio_photo.nml"'
   LEVEL2 'If that is not what you want you have to supply "bio_photo.nml"'
   LEVEL2 'See the bio example on www.gotm.net for a working "bio_photo.nml"'
   return
99 FATAL 'I could not read "bio_photo.nml"'
   stop 'init_bio_photo'

   end subroutine init_bio_photo
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise photo variables
!
! !INTERFACE:
   subroutine init_var_photo
!
! !DESCRIPTION:
! According to the initial water column parameters, the variables of 
! the {\tt photo} module are initialized here. This routine should 
! only be called from inside GOTM. If the {\tt photo} module is 
! linked to an external calling programm (3D model), all variables
! of the model need to initialized from inside the calling program.
!
! !USES:
   IMPLICIT NONE
!

! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Lars Umlauf, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:

   integer                    :: nt,np
   REALTYPE                   :: par

!-----------------------------------------------------------------------
!BOC

   LEVEL2 'init_bio_var'

! fill-in concentration arays and sinking speed
  if (bio_eulerian) then
  
     !  all fields empty
     cc = _ZERO_
     
     !  no sinking
     ws = _ZERO_
     
  else
     
     !  all fields empty
     cc = _ZERO_
     
     !  no sinking
     ws = _ZERO_

     ! initalize particle positions
     call init_par_photo()

     !  initialize particle properties
     nt = 1
     
     do np=1,npar
        
        ! short wave light fraction
        par                  = rad(nlev)*((1.-aa)*exp(par_z(np,nt)/g2)) 
        
        ! local inhibition state      
        par_prop(np,XInd,nt) = _ONE_ - exp(-((max(par,eb)-eb)/eb)**2)        
      
        ! initially fully adapted
        par_prop(np,YInd,nt) =  par_prop(np,XInd,nt)                    

     end do ! particle loop
     
     par_prop(:,ProdInd,nt)  = _ZERO_  
     par_prop(:,CarbInd,nt)  = cinit
   
     
  endif
  

   LEVEL3 'variables initialized'


   return

   end subroutine init_var_photo
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize particle distribution
!
! !INTERFACE:
   subroutine init_par_photo()
!
! !DESCRIPTION:
! Particles are distributed homogeneously over the whole water column. 
! Indices of vertical grid cells are assigned to all particles, and the
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
   end subroutine init_par_photo
!EOC
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update Eulerian variables
!
! !INTERFACE:
   subroutine do_bio_photo_eul()
!
! !DESCRIPTION:
!  Eulerian model is not yet implemented.

! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf,  Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   
   return
   end subroutine do_bio_photo_eul
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update particle properties
!
! !INTERFACE:
   subroutine do_bio_photo_par
!
! !DESCRIPTION: 
! This routine updates the {\tt photo} module according to the updated
! water column properties, i.e.\ after a call to {\tt set\_env\_bio()}. 
! It may be called either from inside GOTM, or from an external 
! calling programm (3D model).
!
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf, Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
   integer                    :: np,nt,i
   REALTYPE                   :: dayhour
   REALTYPE                   :: w,zold,step
   REALTYPE                   :: par,p,pd,pl,x,y,c

!-----------------------------------------------------------------------
!BOC
   
   ! only one type of particles used in this model
   nt = 1

   ! set daytime
   dayhour = secondsofday/secs_pr_hour

   ! vertical migration according to time and inhibition state
   if (photo_mig) then
 
      do np=1,npar
         
         y     = par_prop(np,YInd,nt)   ! get inhibition state
         
         if ((dayhour.gt.6.).and.(dayhour.lt.18.)) then
            
            if ( y.gt.yh )  then
               w = -w_c          ! it's day but it's too bright: move down
            else
               w = w_c           ! it's day and not too bright: move up
            endif
            
         else
            
            w = -w_c             ! it's night: move down
            
         end if
         

!        now move
         step   = dt*w
         zold   = par_z(np,nt)
         par_z(np,nt) = zold + step


!        reflect particle if it jumps below bed or above surface
         if (par_z(np,nt) .lt. zbot) par_z(np,nt) = zbot + ( zbot - par_z(np,nt) )
      
         if (par_z(np,nt) .gt. ztop) par_z(np,nt)=  ztop - par_z(np,nt) 


!        calculate new layer index of particle location
         step=par_z(np,nt) - zold
         if (step.gt.0) then ! search new index above old index
            do i=par_ind(np,nt),nlev
               if (zlev(i) .gt. par_z(np,nt)) EXIT
            end do
         else                ! search new index below old index
            do i=par_ind(np,nt),1,-1
               if (zlev(i-1) .lt. par_z(np,nt)) EXIT
            end do
         end if
         par_ind(np,nt) = i


      end do  ! npar

   end if  ! particle loop


!  compute light inhibition properties, production, and carbon content 

   do np=1,npar
     
      ! short wave light fraction
      par   = rad(nlev)*((_ONE_ - aa)*exp(par_z(np,nt)/g2))    
      
      ! inhibition state
      x     = 1. - exp(-((max(par,eb)-eb)/eb)**2) ! adapted

      y     = par_prop(np,YInd,nt)                ! individual


      ! production
      pd    = pdm*(1.-exp(-par/ed))               ! uninhibited production

      pl    = plm*(1.-exp(-par/el))               ! inhibited production

      y     = y + dt/gamma*(x-y)                  ! update ihibition
      
      p     = pd+y*(pl-pd)                        ! instantaneous production

!     update carbon content
      c     = par_prop(np,CarbInd,nt)     
      c     = c + (p - mu*c)*dt 

!     copy back to property array
      par_prop(np,XInd   ,nt)     = x      
      par_prop(np,YInd   ,nt)     = y      
      par_prop(np,ProdInd,nt)     = p*secs_pr_hour                
      par_prop(np,CarbInd,nt)     = c

   end do ! particle loop

   call do_statistics
   
   return
   end subroutine do_bio_photo_par
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
!  Computes some statistical properties of particles. At the moment
!  only the cell averarages are computed, and cell concentrations are
!  converted to volume concentrations.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf,  Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   integer                    :: ind,np,nt,i

!-----------------------------------------------------------------------
!BOC

   cc = _ZERO_

   nt = 1;   ! only one particle type here

   do np=1,npar
      if (par_act(np,nt)) then

         ind = par_ind(np,nt) 

         ! count particles per grid volume
         cc(1,ind) = cc(1,ind) + _ONE_

         ! add cell properties
         cc(2,ind) = cc(2,ind) + par_prop(np,XInd   ,nt)
         cc(3,ind) = cc(3,ind) + par_prop(np,YInd   ,nt)
         cc(4,ind) = cc(4,ind) + par_prop(np,ProdInd,nt)
         cc(5,ind) = cc(5,ind) + par_prop(np,CarbInd,nt)

      end if

   end do

   ! compute volume averages
   do i=1,nlev
      if (cc(1,i) /= 0) then
         cc(2,i) = cc(2,i)/cc(1,i)
         cc(3,i) = cc(3,i)/cc(1,i)
         cc(4,i) = cc(4,i)/h(i)
         cc(5,i) = cc(5,i)/h(i)
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
   subroutine end_bio_photo
!
! !DESCRIPTION:
!  Nothing done at the moment.
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
   end subroutine end_bio_photo
!EOC

!-----------------------------------------------------------------------

   end module bio_photo

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

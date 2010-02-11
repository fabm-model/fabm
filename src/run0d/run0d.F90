#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 0D shell based on GOTM --- the general framework
!
! !INTERFACE:
   module run0d
!
! !DESCRIPTION:
! TODO
!
! !USES:
   use time
   use bio_0d_base
   use bio_npzd_0d

   IMPLICIT NONE
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_run, time_loop, clean_up

!
! !DEFINED PARAMETERS:
   integer, parameter        :: namlst=10, env_unit=11, out_unit = 12, bio_unit=22
   integer, parameter        :: READ_SUCCESS=1
   integer, parameter        :: END_OF_FILE=-1
   integer, parameter        :: READ_ERROR=-2
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!  private data members initialised via namelists
   character(len=80)         :: title
   REALTYPE                  :: dt
!  station description
   character(len=80)         :: name
   REALTYPE                  :: latitude,longitude

   ! Bio model info
   integer :: bio_model
   logical :: calculate_swr = .true.
   REALTYPE :: cloud
   type (type_model_info) :: modelinfo
   type (type_variable_info), allocatable :: varinfo(:)
   character(len=128)         :: cbuf

   integer, parameter :: npzd_0d_id = 1001
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the model \label{initGOTM}
!
! !INTERFACE:
   subroutine init_run()
!
! !DESCRIPTION:
!  This internal routine triggers the initialization of the model.
!  The first section reads the namelists of {\tt gotmrun.nml} with
!  the user specifications. Then, one by one each of the modules are
!  initialised with help of more specialised routines like
!  {\tt init\_meanflow()} or {\tt init\_turbulence()} defined inside
!  their modules, respectively.
!
!  Note that the KPP-turbulence model requires not only a call to
!  {\tt init\_kpp()} but before also a call to {\tt init\_turbulence()},
!  since there some fields (fluxes, diffusivities, etc) are declared and
!  the turbulence namelist is read.

! !USES:
  IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the gotm module
!
!EOP
!
! !LOCAL VARIABLES:
   character(len=PATH_MAX)   :: env_file,output_file
   integer                   :: i

   namelist /model_setup/ title,dt
   namelist /station/     name,latitude,longitude
   namelist /time/        timefmt,MaxN,start,stop
   namelist /bio/         bio_model,env_file,output_file,calculate_swr,cloud
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_run0d'
   STDERR LINE

!  open the namelist file.
   LEVEL2 'reading model setup namelists..'
   open(namlst,file='run.nml',status='old',action='read',err=90)

   read(namlst,nml=model_setup,err=91)
   read(namlst,nml=station,err=92)
   read(namlst,nml=time,err=93)

   close (namlst)

   LEVEL2 'done.'

!  initialize a few things from  namelists
   timestep   = dt
   !depth0     = depth

!  write information for this run
   LEVEL2 trim(title)
   LEVEL2 'The station ',trim(name),' is situated at (lat,long) ',      &
           latitude,longitude
   LEVEL2 trim(name)

   LEVEL2 'initializing modules....'
   call init_time(MinN,MaxN)

   open(env_unit,file=env_file,action='read', &
        status='old',err=95)
   LEVEL2 'Reading local environment data from:'
   LEVEL3 trim(env_file)
   
   open(out_unit,file=output_file,action='write', &
        status='old',err=96)
   LEVEL2 'Writing results to:'
   LEVEL3 trim(output_file)

   ! Get info on the selected model.
   call init_model_info(modelinfo)
   select case (bio_model)
      case (npzd_0d_id)
         call init_bio_npzd_0d(namlst,'bio_npzd.nml',bio_unit,modelinfo)
   end select
   
   ! Initialize variable information with defaults
   allocate(varinfo(modelinfo%numc))
   do i=1,modelinfo%numc
      call init_variable_info(varinfo(i))
   end do

   ! Get actual variable info from the selected biogeochemical model
   select case (bio_model)
      case (npzd_0d_id)  ! NPZD
         call get_var_info_npzd_0d(modelinfo%numc,varinfo)
   end select

   LEVEL2 'done.'
   STDERR LINE

   return

90 FATAL 'I could not open gotmrun.nml for reading'
   stop 'init_run0d'
91 FATAL 'I could not read the "model_setup" namelist'
   stop 'init_run0d'
92 FATAL 'I could not read the "station" namelist'
   stop 'init_run0d'
93 FATAL 'I could not read the "time" namelist'
   stop 'init_run0d'
94 FATAL 'I could not read the "output" namelist'
   stop 'init_run0d'
95 FATAL 'I could not open ',trim(env_file)
   stop 'init_run0d'
96 FATAL 'I could not open ',trim(output_file)
   stop 'init_run0d'
   end subroutine init_run
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Manage global time--stepping \label{timeLoop}
!
! !INTERFACE:
   subroutine time_loop()
!
! !DESCRIPTION:
! This internal routine is the heart of the code. It contains
! the main time-loop inside of which all routines required
! during the time step are called. The following main processes are
! successively triggered.
! \begin{enumerate}
!  \item The model time is updated and the output is prepared.
!  \item Air-sea interactions (flux, SST) are computed.
!  \item The time step is performed on the mean-flow equations
!        (momentum, temperature).
!  \item Some quantities related to shear and stratification are updated
!        (shear-number, buoyancy frequency, etc).
!  \item Turbulence is updated depending on what turbulence closure
!        model has been specified by the user.
!  \item The results are written to the output files.
! \end{enumerate}
!
! Depending on macros set for the Fortran pre-processor, extra features
! like the effects of sea-grass or sediments are considered in this routine
! (see \sect{sec:extra}).
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the gotm module
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: n
   REALTYPE :: swr
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'time_loop'

   do n=MinN,MaxN

!     prepare time and output
      call update_time(n)
      
      if (calculate_swr) then
         call short_wave_radiation(julianday,secondsofday,longitude,latitude,cloud,swr)
      end if

   end do
   STDERR LINE

   return
   end subroutine time_loop
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read local environment, interpolate in time
!
! !INTERFACE:
   subroutine read_environment(jul,secs,env)
!
! !DESCRIPTION:
!  This routine reads the short wave radiation (in W\,s$^{-2}$) from
!  {\tt swr\_file} and interpolates in time.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: jul,secs
!
! !OUTPUT PARAMETERS:
   type(type_environment),intent(inout) :: env
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,alpha
   REALTYPE, save            :: dt
   integer, save             :: swr_jul1,swr_secs1
   integer, save             :: swr_jul2=0,swr_secs2=0
   REALTYPE, save            :: obs1(1),obs2(1)=0.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   !if (init_saved_vars) then
   !   swr_jul2=0
   !   swr_secs2=0
   !   obs2(1)=0.
   !end if
!  This part initialise and read in new values if necessary.
   if(time_diff(swr_jul2,swr_secs2,jul,secs) .lt. 0) then
      do
         swr_jul1 = swr_jul2
         swr_secs1 = swr_secs2
         obs1 = obs2
         call read_obs(env_unit,yy,mm,dd,hh,min,ss,1,obs2,rc)
         call julian_day(yy,mm,dd,swr_jul2)
         swr_secs2 = hh*3600 + min*60 + ss
         if(time_diff(swr_jul2,swr_secs2,jul,secs) .gt. 0) EXIT
      end do
      dt = time_diff(swr_jul2,swr_secs2,swr_jul1,swr_secs1)
   end if

   call init_environment(env)

!  Do the time interpolation
   t  = time_diff(jul,secs,swr_jul1,swr_secs1)
   alpha = (obs2(1)-obs1(1))/dt
   env%par = obs1(1) + t*alpha

   return
   end subroutine read_environment
   

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_obs
!
! !INTERFACE:
   subroutine read_obs(unit,yy,mm,dd,hh,min,ss,N,obs,ierr)
!
! !DESCRIPTION:
!  This routine will read all non-profile observations.
!  The routine allows for reading more than one scalar variable at a time.
!  The number of data to be read is specified by {\tt N}.
!  Data read-in are returned
!  in the 'obs' array. It is up to the calling routine to assign
!  meaning full variables to the individual elements in {\tt obs}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: unit
   integer, intent(in)                 :: N
!
! !OUTPUT PARAMETERS:
   integer, intent(out)                :: yy,mm,dd,hh,min,ss
   REALTYPE,intent(out)                :: obs(:)
   integer, intent(out)                :: ierr
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See observation module
!
!EOP
!
! !LOCAL VARIABLES:
   character                 :: c1,c2,c3,c4
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
   ierr=0
   read(unit,'(A128)',ERR=100,END=110) cbuf
   read(cbuf,900,ERR=100,END=110) yy,c1,mm,c2,dd,hh,c3,min,c4,ss
   read(cbuf(20:),*,ERR=100,END=110) (obs(i),i=1,N)

   return
100 ierr=READ_ERROR
   return
110 ierr=END_OF_FILE
   return
900 format(i4,a1,i2,a1,i2,1x,i2,a1,i2,a1,i2)
   end subroutine read_obs
!EOC

   
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: The run is over --- now clean up.
!
! !INTERFACE:
   subroutine clean_up()
!
! !DESCRIPTION:
! This function is just a wrapper for the external routine
! {\tt close\_output()} discussed in \sect{sec:output}. All open files
! will be closed after this call.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
!  See log for the gotm module
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'clean_up'
   
   close(env_unit)
   close(out_unit)

   return
   end subroutine clean_up
!EOC

!-----------------------------------------------------------------------

   end module run0d

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

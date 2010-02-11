!$Id: driver_0d.F90,v 1.7 2009-05-10 18:34:51 jorn Exp $
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 0D independent biogeochemical driver based on GOTM
!
! !INTERFACE:
   module driver_0d
!
! !DESCRIPTION:
! TODO
!
! !USES:
   use time
   use bio_types
   use bio_0d_gen

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
   character, parameter      :: separator = char(9)
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
   REALTYPE                  :: latitude,longitude

   ! Bio model info
   integer  :: ode_method, nsave = 1
   integer  :: swr_method = 0
   REALTYPE :: cloud = _ZERO_, depth = _ZERO_, par_fraction = _ONE_, par_background_extinction = _ZERO_
   logical  :: apply_self_shading = .true., add_environment = .false., add_conserved_quantities = .false.
   
   REALTYPE,allocatable                   :: cc(:,:),totals(:)
   type (type_model)                      :: model
   type (type_environment)                :: env
   character(len=128)                     :: cbuf

!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the model
!
! !INTERFACE:
   subroutine init_run()
!
! !DESCRIPTION:
!  This internal routine triggers the initialization of the model.
!  The first section reads the namelists of {\tt run.nml} with
!  the user specifications. Then, one by one each of the modules are
!  initialised.
!
! !USES:
  IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   character(len=PATH_MAX)   :: env_file,output_file
   character(len=64)         :: model_name
   integer                   :: bio_model,i

   namelist /model_setup/ title,bio_model,start,stop,dt,ode_method
   namelist /environment/ env_file,swr_method, &
                          latitude,longitude,cloud,par_fraction, &
                          depth,par_background_extinction,apply_self_shading
   namelist /output/      output_file,nsave,add_environment, &
                          add_conserved_quantities
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_run0d'
   STDERR LINE

   ! Open the namelist file.
   LEVEL2 'reading model setup namelists..'
   open(namlst,file='run.nml',status='old',action='read',err=90)

   ! Read all namelists
   read(namlst,nml=model_setup,err=91)
   read(namlst,nml=environment,err=92)
   read(namlst,nml=output     ,err=93)
   
   ! Close the namelist file.
   close (namlst)

   LEVEL2 'done.'

   ! Configure the time module to use actual start and stop dates.
   timefmt = 2

   ! Transfer the time step to the time module.
   timestep = dt

   ! Write information for this run to the console.
   LEVEL2 trim(title)
   LEVEL2 'The simulated location is situated at (lat,long) ',      &
           latitude,longitude

   LEVEL2 'initializing modules....'

   ! Initialize the time module.
   call init_time(MinN,MaxN)

   ! Open the file with observations of the local environment.
   open(env_unit,file=env_file,action='read', &
        status='old',err=95)
   LEVEL2 'Reading local environment data from:'
   LEVEL3 trim(env_file)
   
   ! Get info on the selected biogeochemical model.
   model_name = get_model_name(bio_model)
   model = init_bio_0d_generic(bio_model,namlst,'bio_'//trim(model_name)//'.nml')
   
   allocate(totals(1:model%info%conserved_quantity_count))

   ! Create state variable vector, using the initial values specified by the model.
   allocate(cc(model%info%state_variable_count,0:1))
   do i=1,model%info%state_variable_count
      cc(i,1) = model%info%variables(i)%initial_value
   end do

   ! Initialize the abiotic environment (this sets all unused variables to a reasonable default)
   call init_environment(env)
   env%z = -depth

   ! Open the output file.
   open(out_unit,file=output_file,action='write', &
        status='replace',err=96)
   LEVEL2 'Writing results to:'
   LEVEL3 trim(output_file)

   ! Write header to the output file.
   write(out_unit,FMT='(''# '',A)') title
   write(out_unit,FMT='(''# '',A)',ADVANCE='NO') 'time'
   if (add_environment) then
      write(out_unit,FMT=100,ADVANCE='NO') separator,'photosynthetically active radiation','W/m2'
      write(out_unit,FMT=100,ADVANCE='NO') separator,'temperature',                        'degrees C'
      write(out_unit,FMT=100,ADVANCE='NO') separator,'salinity',                           'kg/m3'
   end if
   do i=1,model%info%state_variable_count
      write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%info%variables(i)%longname),trim(model%info%variables(i)%unit)
   end do
   if (add_conserved_quantities) then
      do i=1,model%info%conserved_quantity_count
         write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%info%conserved_quantities(i)%longname),trim(model%info%conserved_quantities(i)%unit)
      end do
   end if
   write(out_unit,*)

   LEVEL2 'done.'
   STDERR LINE

   return

90 FATAL 'I could not open run.nml for reading'
   stop 'init_run0d'
91 FATAL 'I could not read the "model_setup" namelist'
   stop 'init_run0d'
92 FATAL 'I could not read the "environment" namelist'
   stop 'init_run0d'
93 FATAL 'I could not read the "output" namelist'
   stop 'init_run0d'
95 FATAL 'I could not open ',trim(env_file)
   stop 'init_run0d'
96 FATAL 'I could not open ',trim(output_file)
   stop 'init_run0d'
100 format (A, A, ' (', A, ')')
   end subroutine init_run
!EOC

!BOP
!
! !IROUTINE: Get the right-hand side of the ODE system.
!
! !INTERFACE:
   subroutine get_rhs(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! TODO
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
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE :: diag(1:model%info%diagnostic_variable_count)
!EOP
!-----------------------------------------------------------------------
!BOC
   call do_bio_0d_generic(model,first,cc(:,1),env,pp(:,:,1),dd(:,:,1),diag)
   
   end subroutine get_rhs
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
! during the time step are called.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i,n
   REALTYPE                  :: extinction
   character(len=19)         :: ts
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'time_loop'

   do n=MinN,MaxN

      ! Update time
      call update_time(n)
      
      ! Update environment
      call read_environment(julianday,secondsofday,env)
      
      ! Calculate photosynthetically active radiation if it is not provided in the input file.
      if (swr_method.eq.0) then
         ! Calculate photosynthetically active radiation from geographic location, time, cloud cover.
         call short_wave_radiation(julianday,secondsofday,longitude,latitude,cloud,env%par)
      end if
      
      ! Apply light attentuation with depth, unless local light is provided in the input file.
      if (swr_method.ne.2) then
         ! Either we calculate surface PAR, or surface PAR is provided.
         ! Calculate the local PAR at the given depth from par fraction, extinction coefficient, and depth.
         extinction = par_background_extinction
         if (apply_self_shading) extinction = extinction + get_bio_extinction_bio_0d_generic(model,cc(:,1))
         env%par = env%par*exp(env%z*extinction)
      end if
      
      ! Multiply by fraction of short-wave radiation that is photosynthetically active.
      env%par = par_fraction*env%par

      ! Integrate one time step
      call ode_solver(ode_method,model%info%state_variable_count,1,dt,cc,get_rhs)
      
      ! Do output
      if (mod(n,nsave).eq.0) then
         call write_time_string(julianday,secondsofday,timestr)
         write (out_unit,FMT='(A)',ADVANCE='NO') timestr
         if (add_environment) then
            write (out_unit,FMT='(A,E14.8E2)',ADVANCE='NO') separator,env%par
            write (out_unit,FMT='(A,E14.8E2)',ADVANCE='NO') separator,env%t
            write (out_unit,FMT='(A,E14.8E2)',ADVANCE='NO') separator,env%s
         end if
         do i=1,model%info%state_variable_count
            write (out_unit,FMT='(A,E14.8E2)',ADVANCE='NO') separator,cc(i,1)
         end do
         if (add_conserved_quantities) then
            totals = get_conserved_quantities_bio_0d_generic(model,cc(:,1))
            do i=1,model%info%conserved_quantity_count
               write (out_unit,FMT='(A,E14.8E2)',ADVANCE='NO') separator,totals(i)
            end do
         end if
         write (out_unit,*)
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
!  This routine reads the local environment from {\tt env\_file} and
!  interpolates in time.
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
   integer,parameter         :: nobs = 3
   integer                   :: yy,mm,dd,hh,min,ss
   REALTYPE                  :: t,curobs(nobs)
   REALTYPE, save            :: dt
   integer, save             :: env_jul1,env_secs1
   integer, save             :: env_jul2=0,env_secs2=0
   REALTYPE, save            :: obs1(nobs),obs2(nobs)=0.
   logical, save             :: endoffile = .false.
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   ! This part reads in new observations if the last read observation lies
   ! before the current time, and the end of the observation file has not yet
   ! been reached.
   if(time_diff(env_jul2,env_secs2,jul,secs) .lt. 0 .and. .not. endoffile) then
      do
         ! Store the previous right-side observation as the new left-side observation.
         env_jul1  = env_jul2
         env_secs1 = env_secs2
         obs1      = obs2
         
         ! Try to read another observation
         call read_obs(env_unit,yy,mm,dd,hh,min,ss,nobs,obs2,rc)
         
         ! Interpret the status code that was returned.
         if (rc==END_OF_FILE) then
            ! Unable to read more observations: set the last observation equal to the first
            ! (the last valid observation read), and stop reading the observation file.
            env_jul2  = env_jul1
            env_secs2 = env_secs1
            obs2      = obs1
            endoffile = .true.
            EXIT
         elseif (rc==READ_ERROR) then
            ! Unknown error occurred: fail.
            stop 'read_environment: error when reading environment from file.'
         end if
         
         ! Calculate the time of the observation we just read.
         call julian_day(yy,mm,dd,env_jul2)
         env_secs2 = hh*3600 + min*60 + ss
         
         ! If the new observation lies beyond the current time, we are done.
         if(time_diff(env_jul2,env_secs2,jul,secs) .gt. 0) EXIT
      end do
      
      if (env_jul1.eq.0) then
         ! The time of the very first observation already lies beyond current time.
         ! Set the left-side observation equal to the right-side one that we just read.
         env_jul1  = env_jul2
         env_secs1 = env_secs2
         obs1      = obs2
      end if
      
      ! Calculate the difference in time between the left- and right-side observation.
      dt = time_diff(env_jul2,env_secs2,env_jul1,env_secs1)
   end if

   if (dt.eq.0) then
      ! The time of both observations is identical: we could not get one observation
      ! before and one observation after the current time. Do not interpolate and just use
      ! the only [nearest] observation as-is.
      curobs = obs1
   else
      ! Interpolate in time
      t  = time_diff(jul,secs,env_jul1,env_secs1)
      curobs = obs1 + t*(obs2-obs1)/dt
   end if

   ! Store environment properties.
   env%par = curobs(1)
   env%t   = curobs(2)
   env%s   = curobs(3)

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
! Close all open files.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
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

   end module driver_0d

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

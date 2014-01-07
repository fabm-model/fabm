#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 0D independent driver for the Framework for Aquatic Biogeochemical Models (FABM)
!
! !INTERFACE:
   module fabm_0d
!
! !DESCRIPTION:
! TODO
!
! !USES:
   use time
   use input
   use fabm
   use fabm_types, only:rk

   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_run, time_loop, clean_up

!
! !DEFINED PARAMETERS:
   integer, parameter        :: namlst=10, out_unit = 12, bio_unit=22
   integer, parameter        :: READ_SUCCESS=1
   integer, parameter        :: END_OF_FILE=-1
   integer, parameter        :: READ_ERROR=-2
   integer, parameter        :: CENTER=0,SURFACE=1,BOTTOM=2
   character, parameter      :: separator = char(9)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!  private data members initialised via namelists
   character(len=80)         :: title
   real(rk)                  :: dt
!  station description
   real(rk)                  :: latitude,longitude

   ! Bio model info
   integer  :: ode_method
   logical  :: repair_state
   integer(timestepkind)::nsave
   integer  :: swr_method
   real(rk) :: cloud
   real(rk) :: par_fraction
   real(rk) :: par_background_extinction
   logical  :: apply_self_shading
   logical  :: add_environment
   logical  :: add_conserved_quantities
   logical  :: add_diagnostic_variables

   ! Environment
   real(rk),target :: temp,salt,par,current_depth,dens,wind_sf,taub,decimal_yearday
   real(rk)        :: par_sf,par_bt,par_ct,column_depth

   real(rk),allocatable      :: cc(:,:),totals(:)
   type (type_model),pointer :: model
   character(len=128)        :: cbuf

   interface
      function short_wave_radiation(jul,secs,dlon,dlat,cloud,bio_albedo) result(swr)
         import
         integer, intent(in)                 :: jul,secs
         real(rk), intent(in)                :: dlon,dlat
         real(rk), intent(in)                :: cloud
         real(rk), intent(in)                :: bio_albedo
         real(rk)                            :: swr
      end function short_wave_radiation
   end interface
!EOP
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   character(len=PATH_MAX)   :: env_file,output_file
   integer                   :: i
   real(rk)                  :: depth
   real(rk),parameter        :: invalid_latitude = -100._rk,invalid_longitude = -400.0_rk

   namelist /model_setup/ title,start,stop,dt,ode_method,repair_state
   namelist /environment/ env_file,swr_method, &
                          latitude,longitude,cloud,par_fraction, &
                          depth,par_background_extinction,apply_self_shading
   namelist /output/      output_file,nsave,add_environment, &
                          add_diagnostic_variables, add_conserved_quantities
!
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'init_run'
   STDERR LINE

   ! Open the namelist file.
   LEVEL2 'reading model setup namelists..'
   open(namlst,file='run.nml',status='old',action='read',err=90)

   ! Initialize environment
   temp = 0.0_rk
   salt = 0.0_rk
   par = 0.0_rk
   dens = 0.0_rk
   wind_sf = 0.0_rk
   par_sf = 0.0_rk
   par_bt = 0.0_rk
   par_ct = 0.0_rk
   decimal_yearday = 0.0_rk
   taub  = 0.0_rk

   ! Read all namelists
   title = ''
   start = ''
   stop = ''
   dt = 0.0_rk
   ode_method = 1
   repair_state = .false.
   read(namlst,nml=model_setup,err=91)
   
   ! Read environment namelist
   env_file = ''
   swr_method = 0
   latitude = invalid_latitude
   longitude = invalid_longitude
   cloud = 0.0_rk
   par_fraction = 1.0_rk
   depth = -1.0_rk
   par_background_extinction = 0.0_rk
   apply_self_shading = .true.
   read(namlst,nml=environment,err=92)

   ! Read output namelist
   output_file = ''
   nsave = 1
   add_environment = .false.
   add_conserved_quantities = .false.
   add_diagnostic_variables=.false.
   read(namlst,nml=output     ,err=93)

   ! Close the namelist file.
   close (namlst)

   if (start=='') then
      FATAL 'run.nml: start time "start" must be set in "model_setup" namelist.'
      stop 'init_run'
   end if

   if (stop=='') then
      FATAL 'run.nml: stop time "stop" must be set in "model_setup" namelist.'
      stop 'init_run'
   end if

   if (dt<=0.0_rk) then
      FATAL 'run.nml: time step "dt" must be set to a positive value in "model_setup" namelist.'
      stop 'init_run'
   end if

   if (env_file=='') then
      FATAL 'run.nml: "env_file" must be set to a valid file path in "environment" namelist.'
      stop 'init_run'
   end if

   if (latitude/=invalid_latitude.and.(latitude<-90._rk.or.latitude>90._rk)) then
      FATAL 'run.nml: latitude must lie between -90 and 90.'
      stop 'init_run'
   end if
   if (longitude/=invalid_longitude.and.(longitude<-360._rk.or.longitude>360._rk)) then
      FATAL 'run.nml: longitude must lie between -360 and 360.'
      stop 'init_run'
   end if

   ! Make sure depth has been provided.
   if (depth<=0.0_rk) then
      FATAL 'run.nml: a positive value for "depth" must be provided in "environment" namelist.'
      stop 'init_run'
   end if
   column_depth = depth ! For now we have a single depth value only. Use that for both column depth and evaluation depth.
   call update_depth(CENTER)
   
   ! If longitude and latitude are used, make sure they have been provided and are valid.
   if (swr_method==0) then
      if (latitude==invalid_latitude) then
         FATAL 'run.nml: a valid value for "latitude" must be provided in "environment" if "swr_method" is 0.'
         stop 'init_run'
      end if
      if (longitude==invalid_longitude) then
         FATAL 'run.nml: a valid value for "longitude" must be provided in "environment" if "swr_method" is 0.'
         stop 'init_run'
      end if
   end if
   
   if (output_file=='') then
      FATAL 'run.nml: "output_file" must be set to a valid file path in "output" namelist.'
      stop 'init_run'
   end if

   LEVEL2 'done.'

   ! Configure the time module to use actual start and stop dates.
   timefmt = 2

   ! Transfer the time step to the time module.
   timestep = dt

   ! Write information for this run to the console.
   LEVEL2 'Simulation: '//trim(title)
   select case (swr_method)
      case (0)
         LEVEL2 'Surface photosynthetically active radiation will be calculated from time,'
         LEVEL2 'cloud cover, and the simulated location at (lat,long)'
         LEVEL2 latitude,longitude
         LEVEL2 'Local PAR will be calculated from the surface value,'
         LEVEL2 'depth, and light extinction coefficient.'
      case (1)
         LEVEL2 'Surface photosynthetically active radiation (PAR) is provided as input.'
         LEVEL2 'Local PAR will be calculated from the surface value,'
         LEVEL2 'depth, and light extinction coefficient.'
      case (2)
         LEVEL2 'Local photosynthetically active radiation is provided as input.'
   end select

   LEVEL2 'initializing modules....'

   ! Initialize the time module.
   call init_time(MinN,MaxN)

   ! Open the file with observations of the local environment.
   LEVEL2 'Reading local environment data from:'
   LEVEL3 trim(env_file)
   call init_input()
   call register_input_0d(env_file,1,par_sf)
   call register_input_0d(env_file,2,temp)
   call register_input_0d(env_file,3,salt)

   ! Build FABM model tree.
   model => fabm_create_model_from_file(namlst)

   ! Send information on spatial domain to FABM (this also allocates memory for diagnostics)
   call fabm_set_domain(model,dt)

   ! Allocate space for totals of conserved quantities.
   allocate(totals(1:size(model%info%conserved_quantities)))

   ! Create state variable vector, using the initial values specified by the model,
   ! and link state data to FABM.
   allocate(cc(size(model%info%state_variables)+size(model%info%state_variables_ben),0:1))
   do i=1,size(model%info%state_variables)
      cc(i,1) = model%info%state_variables(i)%initial_value
      call fabm_link_bulk_state_data(model,i,cc(i,1))
   end do

   ! Create benthos variable vector, using the initial values specified by the model,
   ! and link state data to FABM.
   do i=1,size(model%info%state_variables_ben)
      cc(size(model%info%state_variables)+i,1) = model%info%state_variables_ben(i)%initial_value
      call fabm_link_bottom_state_data(model,i,cc(size(model%info%state_variables)+i,1))
   end do

   ! Link environmental data to FABM
   call fabm_link_bulk_data(model,standard_variables%temperature,temp)
   call fabm_link_bulk_data(model,standard_variables%practical_salinity,salt)
   call fabm_link_bulk_data(model,standard_variables%downwelling_photosynthetic_radiative_flux,par)
   call fabm_link_bulk_data(model,standard_variables%pressure,current_depth)
   call fabm_link_bulk_data(model,standard_variables%density,dens)
   call fabm_link_bulk_data(model,standard_variables%cell_thickness,column_depth)
   call fabm_link_bulk_data(model,standard_variables%depth,current_depth)
   call fabm_link_horizontal_data(model,standard_variables%wind_speed,wind_sf)
   call fabm_link_horizontal_data(model,standard_variables%surface_downwelling_photosynthetic_radiative_flux,par_sf)
   call fabm_link_horizontal_data(model,standard_variables%bottom_stress,taub)
   call fabm_link_horizontal_data(model,standard_variables%cloud_area_fraction,cloud)
   call fabm_link_horizontal_data(model,standard_variables%bottom_depth,column_depth)
   call fabm_link_horizontal_data(model,standard_variables%bottom_depth_below_geoid,column_depth)
   if (latitude /=invalid_latitude ) call fabm_link_horizontal_data(model,standard_variables%latitude,latitude)
   if (longitude/=invalid_longitude) call fabm_link_horizontal_data(model,standard_variables%longitude,longitude)
   call fabm_link_scalar_data(model,standard_variables%number_of_days_since_start_of_the_year,decimal_yearday)

   call fabm_check_ready(model)

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
   do i=1,size(model%info%state_variables)
      write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%info%state_variables(i)%long_name),trim(model%info%state_variables(i)%units)
   end do
   do i=1,size(model%info%state_variables_ben)
      write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%info%state_variables_ben(i)%long_name),trim(model%info%state_variables_ben(i)%units)
   end do
   if (add_diagnostic_variables) then
      do i=1,size(model%info%diagnostic_variables)
         write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%info%diagnostic_variables(i)%long_name),trim(model%info%diagnostic_variables(i)%units)
      end do
      do i=1,size(model%info%diagnostic_variables_hz)
         write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%info%diagnostic_variables_hz(i)%long_name),trim(model%info%diagnostic_variables_hz(i)%units)
      end do
   end if
   if (add_conserved_quantities) then
      do i=1,size(model%info%conserved_quantities)
         write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%info%conserved_quantities(i)%long_name),trim(model%info%conserved_quantities(i)%units)
      end do
   end if
   write(out_unit,*)

   LEVEL2 'done.'
   STDERR LINE

   return

90 FATAL 'I could not open run.nml for reading'
   stop 'init_run'
91 FATAL 'I could not read the "model_setup" namelist'
   stop 'init_run'
92 FATAL 'I could not read the "environment" namelist'
   stop 'init_run'
93 FATAL 'I could not read the "output" namelist'
   stop 'init_run'
95 FATAL 'I could not open ',trim(env_file)
   stop 'init_run'
96 FATAL 'I could not open ',trim(output_file)
   stop 'init_run'
97 FATAL 'I could not open fabm.nml for reading'
   stop 'init_run'
100 format (A, A, ' (', A, ')')
   end subroutine init_run
!EOC

   subroutine update_depth(location)
      integer, intent(in) :: location
      
      select case (location)
         case (SURFACE)
            current_depth = 0.0_rk
            par = par_sf
         case (BOTTOM)
            current_depth = column_depth            
            par = par_bt
         case (CENTER)
            current_depth = 0.5_rk*column_depth
            par = par_ct
      end select
   end subroutine update_depth

!BOP
!
! !IROUTINE: Get the right-hand side of the ODE system.
!
! !INTERFACE:
   subroutine get_ppdd(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! TODO
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: numc,nlev
   real(rk), intent(in)                 :: cc(1:numc,0:nlev)
!
! !OUTPUT PARAMETERS:
   real(rk), intent(out)                :: pp(1:numc,1:numc,0:nlev)
   real(rk), intent(out)                :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL PARAMETERS:
   integer                              :: n
!EOP

!-----------------------------------------------------------------------
!BOC
   ! Initialize production/destruction matrices to zero (entries will be incremented by FABM)
   pp(:,:,1) = 0.0_rk
   dd(:,:,1) = 0.0_rk

   ! Send current state to FABM
   ! (this may differ from the global state cc if using a multi-step integration scheme such as Runge-Kutta)
   do n=1,size(model%info%state_variables)
      call fabm_link_bulk_state_data(model,n,cc(n,1))
   end do
   do n=1,size(model%info%state_variables_ben)
      call fabm_link_bottom_state_data(model,n,cc(size(model%info%state_variables)+n,1))
   end do

   ! Shortcut to the number of pelagic state variables.
   n = size(model%info%state_variables)

   ! Calculate temporal derivatives due to benthic processes.
   call update_depth(BOTTOM)
   call fabm_do_benthos(model,pp(:,:,1),dd(:,:,1),n)

   ! For pelagic variables: translate bottom flux to into change in concentration
   pp(1:n,:,1) = pp(1:n,:,1)/column_depth
   dd(1:n,:,1) = dd(1:n,:,1)/column_depth

   ! For pelagic variables: surface and bottom flux (rate per surface area) to concentration (rate per volume)
   call update_depth(CENTER)
   call fabm_do(model,pp(:,:,1),dd(:,:,1))

   end subroutine get_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the right-hand side of the ODE system.
!
! !INTERFACE:
   subroutine get_rhs(first,numc,nlev,cc,rhs)
!
! !DESCRIPTION:
! TODO
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: numc,nlev
   real(rk), intent(in)                 :: cc(1:numc,0:nlev)
!
! !OUTPUT PARAMETERS:
   real(rk), intent(out)                :: rhs(1:numc,0:nlev)
!
! !LOCAL PARAMETERS:
   integer                              :: n
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Initialize derivatives to zero (entries will be incremented by FABM)
   rhs(:,1) = 0.0_rk

   ! Send current state to FABM
   ! (this may differ from the global state cc if using a multi-step integration scheme such as Runge-Kutta)
   do n=1,size(model%info%state_variables)
      call fabm_link_bulk_state_data(model,n,cc(n,1))
   end do
   do n=1,size(model%info%state_variables_ben)
      call fabm_link_bottom_state_data(model,n,cc(size(model%info%state_variables)+n,1))
   end do

   ! Shortcut to the number of pelagic state variables.
   n = size(model%info%state_variables)

   ! Calculate temporal derivatives due to surface exchange.
   call update_depth(SURFACE)
   call fabm_get_surface_exchange(model,rhs(1:n,1))

   ! Calculate temporal derivatives due to benthic processes.
   call update_depth(BOTTOM)
   call fabm_do_benthos(model,rhs(1:n,1),rhs(n+1:,1))

   ! For pelagic variables: surface and bottom flux (rate per surface area) to concentration (rate per volume)
   rhs(1:n,1) = rhs(1:n,1)/column_depth

   ! Add change in pelagic variables.
   call update_depth(CENTER)
   call fabm_do(model,rhs(:,1))

   end subroutine get_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check the current values of all state variables
!
! !INTERFACE:
   subroutine do_repair_state(location)
!
! !DESCRIPTION:
!  Checks the current values of all state variables and repairs these
!  if allowed and possible. If the state is invalid and repair is not
!  allowed, the model is brought down.
!
! !INPUT PARAMETERS:
   character(len=*),intent(in) :: location
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
   logical :: valid
!
!-----------------------------------------------------------------------
!BOC

   call fabm_check_state(model,repair_state,valid)
   if (.not. (valid .or. repair_state)) then
      FATAL 'State variable values are invalid and repair is not allowed.'
      FATAL location
      FATAL 'note that repair_state() should be used with caution.'
      FATAL 'try and decrease dt first - and see if that helps.'
      stop 'od_fabm::do_repair_state'
   end if

   end subroutine do_repair_state
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                   :: i
   integer(timestepkind)     :: n
   real(rk)                  :: extinction,bio_albedo
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'time_loop'

   do n=MinN,MaxN

      ! Update time
      call update_time(n)

      decimal_yearday = yearday-1 + dble(secondsofday)/86400.

      ! Update environment
      call do_input(julianday,secondsofday)

      ! Calculate photosynthetically active radiation if it is not provided in the input file.
      if (swr_method==0) then
         ! Calculate photosynthetically active radiation from geographic location, time, cloud cover.
         call fabm_get_albedo(model,bio_albedo)
         par_sf = short_wave_radiation(julianday,secondsofday,longitude,latitude,cloud,bio_albedo)
      end if

      ! Multiply by fraction of short-wave radiation that is photosynthetically active.
      par_sf = par_fraction*par_sf

      ! Apply light attentuation with depth, unless local light is provided in the input file.
      if (swr_method/=2) then
         ! Either we calculate surface PAR, or surface PAR is provided.
         ! Calculate the local PAR at the given depth from par fraction, extinction coefficient, and depth.
         extinction = 0.0_rk
         if (apply_self_shading) call fabm_get_light_extinction(model,extinction)
         extinction = extinction + par_background_extinction
         par_ct = par_sf*exp(-0.5_rk*column_depth*extinction)
         par_bt = par_sf*exp(-column_depth*extinction)
      else
         par_ct = par_sf
         par_bt = par_sf
      end if
      call update_depth(CENTER)

     ! Repair state before calling FABM
     call do_repair_state('0d::time_loop(), before ode_solver()')

     call fabm_update_time(model,real(n,rk))

      ! Integrate one time step
      call ode_solver(ode_method,size(model%info%state_variables)+size(model%info%state_variables_ben),1,dt,cc,get_rhs,get_ppdd)

      ! ODE solver may have redirected the current state with to an array with intermediate values.
      ! Reset to global array.
      do i=1,size(model%info%state_variables)
         call fabm_link_bulk_state_data(model,i,cc(i,1))
      end do
      do i=1,size(model%info%state_variables_ben)
         call fabm_link_bottom_state_data(model,i,cc(size(model%info%state_variables)+i,1))
      end do

     call do_repair_state('0d::time_loop(), after ode_solver()')

      ! Do output
      if (mod(n,nsave)==0) then
         call write_time_string(julianday,secondsofday,timestr)
         write (out_unit,FMT='(A)',ADVANCE='NO') timestr
         if (add_environment) then
            write (out_unit,FMT='(A,E15.8E3)',ADVANCE='NO') separator,par
            write (out_unit,FMT='(A,E15.8E3)',ADVANCE='NO') separator,temp
            write (out_unit,FMT='(A,E15.8E3)',ADVANCE='NO') separator,salt
         end if
         do i=1,(size(model%info%state_variables)+size(model%info%state_variables_ben))
            write (out_unit,FMT='(A,E15.8E3)',ADVANCE='NO') separator,cc(i,1)
         end do
         if (add_diagnostic_variables) then
            do i=1,size(model%info%diagnostic_variables)
               write (out_unit,FMT='(A,E15.8E3)',ADVANCE='NO') separator,fabm_get_bulk_diagnostic_data(model,i)
            end do
            do i=1,size(model%info%diagnostic_variables_hz)
               write (out_unit,FMT='(A,E15.8E3)',ADVANCE='NO') separator,fabm_get_horizontal_diagnostic_data(model,i)
            end do
         end if
         if (add_conserved_quantities) then
            call fabm_get_conserved_quantities(model,totals)
            do i=1,size(model%info%conserved_quantities)
               write (out_unit,FMT='(A,E15.8E3)',ADVANCE='NO') separator,totals(i)
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
! !IROUTINE: The run is over --- now clean up.
!
! !INTERFACE:
   subroutine clean_up()
!
! !DESCRIPTION:
! Close all open files.
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'clean_up'

   call close_input()
   close(out_unit)

   end subroutine clean_up
!EOC

!-----------------------------------------------------------------------

   end module fabm_0d

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

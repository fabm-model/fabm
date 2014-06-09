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

   use fabm_driver
   use fabm_types
   use fabm_expressions
   use fabm
   use fabm_config
#ifdef NETCDF4
   use netcdf
#endif

   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_run, time_loop, clean_up
!
! !DEFINED PARAMETERS:
   integer, parameter        :: namlst=10, out_unit = 12, bio_unit=22, yaml_unit=23
   integer, parameter        :: CENTER=0,SURFACE=1,BOTTOM=2
   character, parameter      :: separator = char(9)
   integer, parameter        :: ASCII_FMT=1
   integer, parameter        :: NETCDF_FMT=2
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!  private data members initialised via namelists
   character(len=80)         :: title
   real(rk)                  :: dt
!  station description
   real(rk)                  :: latitude,longitude
   integer                   :: output_format
   integer                   :: ncid=-1
#ifdef NETCDF4
   integer                   :: setn
   integer                   :: time_id
   integer                   :: par_id,temp_id,salt_id
   integer, allocatable, dimension(:)  :: statevar_ids,diagnostic_ids,conserved_ids,conserved_change_ids
#endif

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
   real(rk),target :: temp,salt,par,current_depth,dens,decimal_yearday
   real(rk)        :: par_sf,par_bt,par_ct,column_depth

   real(rk),allocatable      :: cc(:,:),totals(:),totals0(:),expression_data(:)
   type (type_model),pointer :: model
   character(len=128)        :: cbuf

   type type_input_data
      character(len=attribute_length) :: variable_name = ''
      real(rk)                        :: value         = 0.0_rk
      type (type_input_data),pointer  :: next          => null()
   end type
   type (type_input_data), pointer, save :: first_input_data => null()

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
!
! !LOCAL VARIABLES:
   character(len=PATH_MAX)   :: env_file,output_file
   integer                   :: i
   real(rk)                  :: depth
   real(rk),parameter        :: invalid_latitude = -100._rk,invalid_longitude = -400.0_rk
   logical                   :: file_exists
   real(rk),allocatable      :: rhs(:,:)

   namelist /model_setup/ title,start,stop,dt,ode_method,repair_state
   namelist /environment/ env_file,swr_method, &
                          latitude,longitude,cloud,par_fraction, &
                          depth,par_background_extinction,apply_self_shading
   namelist /output/      output_file,output_format,nsave,add_environment, &
                          add_diagnostic_variables, add_conserved_quantities
!EOP
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
   dens = 1027.0_rk
   par_sf = 0.0_rk
   par_bt = 0.0_rk
   par_ct = 0.0_rk
   decimal_yearday = 0.0_rk

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
   output_format=ASCII_FMT
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
#if 0
   if (output_format/=ASCII_FMT) then
      STDERR 'run.nml: "output_format" NetCDF is not fully implemented yet'
      STDERR 'shifting to ASCII output.'
      output_format=ASCII_FMT
   end if
#endif


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

   ! Build FABM model tree. Use fabm.yaml if available, otherwise fall back to fabm.nml.
   inquire(file='fabm.yaml',exist=file_exists)
   if (file_exists) then
      ! From YAML file fabm.yaml
      allocate(model)
      call fabm_create_model_from_yaml_file(model)
   else
      ! From namelists in fabm.nml
      model => fabm_create_model_from_file(namlst)
   end if

   ! Send information on spatial domain to FABM (this also allocates memory for diagnostics)
   call fabm_set_domain(model,dt)

   ! Allocate space for totals of conserved quantities.
   allocate(totals0(1:size(model%conserved_quantities)))  ! at initial time
   allocate(totals(1:size(model%conserved_quantities)))   ! at current time

   ! Create state variable vector, using the initial values specified by the model,
   ! and link state data to FABM.
   allocate(cc(size(model%state_variables)+size(model%bottom_state_variables),0:1))
   do i=1,size(model%state_variables)
      cc(i,1) = model%state_variables(i)%initial_value
      call fabm_link_bulk_state_data(model,i,cc(i,1))
   end do

   ! Create benthos variable vector, using the initial values specified by the model,
   ! and link state data to FABM.
   do i=1,size(model%bottom_state_variables)
      cc(size(model%state_variables)+i,1) = model%bottom_state_variables(i)%initial_value
      call fabm_link_bottom_state_data(model,i,cc(size(model%state_variables)+i,1))
   end do

   ! Link environmental data to FABM
   call fabm_link_bulk_data(model,standard_variables%temperature,temp)
   call fabm_link_bulk_data(model,standard_variables%practical_salinity,salt)
   call fabm_link_bulk_data(model,standard_variables%downwelling_photosynthetic_radiative_flux,par)
   call fabm_link_bulk_data(model,standard_variables%pressure,current_depth)
   call fabm_link_bulk_data(model,standard_variables%density,dens)
   call fabm_link_bulk_data(model,standard_variables%cell_thickness,column_depth)
   call fabm_link_bulk_data(model,standard_variables%depth,current_depth)
   call fabm_link_horizontal_data(model,standard_variables%surface_downwelling_photosynthetic_radiative_flux,par_sf)
   call fabm_link_horizontal_data(model,standard_variables%cloud_area_fraction,cloud)
   call fabm_link_horizontal_data(model,standard_variables%bottom_depth,column_depth)
   call fabm_link_horizontal_data(model,standard_variables%bottom_depth_below_geoid,column_depth)
   if (latitude /=invalid_latitude ) call fabm_link_horizontal_data(model,standard_variables%latitude,latitude)
   if (longitude/=invalid_longitude) call fabm_link_horizontal_data(model,standard_variables%longitude,longitude)
   call fabm_link_scalar_data(model,standard_variables%number_of_days_since_start_of_the_year,decimal_yearday)

   ! Read forcing data specified in input.yaml.
   call init_yaml()

   ! Allocate memory for the value of any requested vertical integrals/averages of FABM variables.
   call check_fabm_expressions()

   call fabm_check_ready(model)

   ! Allow the model to compute all diagnostics, so output for initial time contains sensible values.
   call update_time(0_timestepkind)
   decimal_yearday = yearday-1 + dble(secondsofday)/86400
   call do_input(julianday,secondsofday)
   allocate(rhs(size(cc,1),0:1))
   call get_rhs(.true.,size(cc,1),1,cc,rhs)

   call init_output(output_format,output_file,start)
   call do_output(output_format,0_timestepkind)

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
97 FATAL 'I could not open fabm.nml for reading'
   stop 'init_run'
   end subroutine init_run
!EOC

   subroutine init_yaml()
      use fabm_config_types
      use fabm_yaml,yaml_parse=>parse,yaml_error_length=>error_length

      logical                            :: exists
      character(len=yaml_error_length)   :: yaml_error
      class (type_node),         pointer :: node

      ! Determine whether input.yaml exists.
      inquire(file='input.yaml',exist=exists)
      if (.not.exists) return

      ! Parse YAML.
      node => yaml_parse('input.yaml',yaml_unit,yaml_error)
      if (yaml_error/='') call driver%fatal_error('init_yaml',trim(yaml_error))

      ! Process root-level dictionary.
      select type (node)
         class is (type_dictionary)
            ! Process "input" section.
            call init_yaml_input(node)
         class default
            call fatal_error('init_yaml','input.yaml must contain a dictionary with (variable name : information) pairs, not a single value.')
      end select
   end subroutine init_yaml

   subroutine init_yaml_input(mapping)
      use fabm_config_types

      class (type_dictionary),intent(in)  :: mapping

      character(len=64)                  :: variable_name
      type (type_key_value_pair),pointer :: pair

      pair => mapping%first
      do while (associated(pair))
         variable_name = trim(pair%key)
         select type (dict=>pair%value)
            class is (type_dictionary)
               call parse_input_variable(variable_name,dict)
            class default
               call fatal_error('init_input','Contents of '//trim(dict%path)//' must be a dictionary, not a single value.')
         end select
         pair => pair%next
      end do
   end subroutine init_yaml_input

   subroutine parse_input_variable(variable_name,mapping)
      use fabm_config_types

      character(len=*),      intent(in) :: variable_name
      type (type_dictionary),intent(in) :: mapping

      type (type_error), pointer        :: config_error
      class (type_node),  pointer :: node
      class (type_scalar),pointer :: scalar
      real(rk)            :: relaxation_time
      character(len=1024) :: path
      integer             :: column
      logical             :: is_state_variable
      type (type_input_data),pointer :: input_data
      type (type_key_value_pair),pointer :: pair

      type (type_bulk_variable_id)       :: bulk_id
      type (type_horizontal_variable_id) :: horizontal_id
      type (type_scalar_variable_id)     :: scalar_id

      ! Try to locate the forced variable among bulk, horizontal, and global variables in the active biogeochemical models.
      is_state_variable = .false.
      bulk_id = fabm_get_bulk_variable_id(model,variable_name)
      if (fabm_is_variable_used(bulk_id)) then
         is_state_variable = bulk_id%state_index/=-1
      else
         horizontal_id = fabm_get_horizontal_variable_id(model,variable_name)
         if (fabm_is_variable_used(horizontal_id)) then
            is_state_variable = horizontal_id%state_index/=-1
         else
            scalar_id = fabm_get_scalar_variable_id(model,variable_name)
            if (.not.fabm_is_variable_used(scalar_id)) &
               call fatal_error('parse_input_variable','input.yaml: Variable "'//trim(variable_name)//'" is not present in any biogeochemical  model.')
         end if
      end if

      ! Create an object to hold information on this input variable.      
      if (associated(first_input_data)) then
         input_data => first_input_data
         do while (associated(input_data%next))
            input_data => input_data%next
         end do
         allocate(input_data%next)
         input_data => input_data%next
      else
         allocate(first_input_data)
         input_data => first_input_data
      end if
      input_data%variable_name = variable_name

      scalar => mapping%get_scalar('constant_value',required=.false.,error=config_error)
      if (associated(scalar)) then
         ! Input variable is set to a constant value.
         input_data%value = mapping%get_real('constant_value',error=config_error)
         if (associated(config_error)) call fatal_error('parse_input_variable',config_error%message)

         node => mapping%get('file')
         if (associated(node)) call fatal_error('parse_input_variable','input.yaml, variable "'//trim(variable_name)//'": keys "constant_value" and "file" cannot both be present.')
         node => mapping%get('column')
         if (associated(node)) call fatal_error('parse_input_variable','input.yaml, variable "'//trim(variable_name)//'": keys "constant_value" and "column" cannot both be present.')
      else
         ! Input variable is set to a time-varying value. Obtain path and column number.
         path = mapping%get_string('file',error=config_error)
         if (associated(config_error)) call fatal_error('parse_input_variable',config_error%message)
         column = mapping%get_integer('column',default=1,error=config_error)
         if (associated(config_error)) call fatal_error('parse_input_variable',config_error%message)
         call register_input_0d(path,column,input_data%value)
      end if

      if (is_state_variable) then
         relaxation_time = mapping%get_real('relaxation_time',default=1.e15_rk,error=config_error)
         if (associated(config_error)) call fatal_error('parse_input_variable',config_error%message)
      else
         node => mapping%get('relaxation_time')
         if (associated(node)) call fatal_error('parse_input_variable','input.yaml, variable "'//trim(variable_name)//'": key "relaxation_time" is not supported because "'//trim(variable_name)//'" is not a state variable.')
      end if

      ! Warn about uninterpreted keys.
      pair => mapping%first
      do while (associated(pair))
         if (.not.pair%accessed) call fatal_error('parse_input_variable','input.yaml: Unrecognized option "'//trim(pair%key)//'" found for variable "'//trim(variable_name)//'".')
         pair => pair%next
      end do

      ! Link forced data to target variable.
      if (fabm_is_variable_used(bulk_id)) then
         call fabm_link_bulk_data(model,bulk_id,input_data%value)
      elseif (fabm_is_variable_used(horizontal_id)) then
         call fabm_link_horizontal_data(model,horizontal_id,input_data%value)
      else
         call fabm_link_scalar_data(model,scalar_id,input_data%value)
      end if

   end subroutine parse_input_variable

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

      decimal_yearday = yearday-1 + dble(secondsofday)/86400

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
      call ode_solver(ode_method,size(model%state_variables)+size(model%bottom_state_variables),1,dt,cc,get_rhs,get_ppdd)

      ! ODE solver may have redirected the current state with to an array with intermediate values.
      ! Reset to global array.
      do i=1,size(model%state_variables)
         call fabm_link_bulk_state_data(model,i,cc(i,1))
      end do
      do i=1,size(model%bottom_state_variables)
         call fabm_link_bottom_state_data(model,i,cc(size(model%state_variables)+i,1))
      end do

      call do_repair_state('0d::time_loop(), after ode_solver()')

      ! Do output
      if (mod(n,nsave)==0) then
         call do_output(output_format,n)
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
! !LOCAL PARAMETERS:
   integer                              :: iret
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'clean_up'

   call close_input()
   if (ncid .ne. -1) then
#ifdef NETCDF4
      iret = nf90_close(ncid)
      call check_err(iret)
#endif
   else
      close(out_unit)
   end if

   end subroutine clean_up
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the right-hand side of the ODE system.
!
! !INTERFACE:
   subroutine update_depth(location)
!
! !DESCRIPTION:
! TODO
!
! !INPUT PARAMETERS:
   integer, intent(in)                  :: location
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
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
   do n=1,size(model%state_variables)
      call fabm_link_bulk_state_data(model,n,cc(n,1))
   end do
   do n=1,size(model%bottom_state_variables)
      call fabm_link_bottom_state_data(model,n,cc(size(model%state_variables)+n,1))
   end do

   call update_fabm_expressions()

   ! Shortcut to the number of pelagic state variables.
   n = size(model%state_variables)

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
   do n=1,size(model%state_variables)
      call fabm_link_bulk_state_data(model,n,cc(n,1))
   end do
   do n=1,size(model%bottom_state_variables)
      call fabm_link_bottom_state_data(model,n,cc(size(model%state_variables)+n,1))
   end do

   call update_fabm_expressions()

   ! Shortcut to the number of pelagic state variables.
   n = size(model%state_variables)

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
! !IROUTINE: prepare for output
!
! !INTERFACE:
   subroutine init_output(output_format,output_file,start)
!
! !DESCRIPTION:
! TODO
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: output_format
   character(len=*), intent(in)        :: output_file,start
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL PARAMETERS:
   integer         :: i,iret
#ifdef NETCDF4
   integer         :: lon_dim,lat_dim,time_dim
   integer         :: lon_id,lat_id
   integer         :: dims(3)
   character(len=128) :: ncdf_time_str
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (output_format)
      case(ASCII_FMT)
         open(out_unit,file=trim(output_file),action='write', &
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
         do i=1,size(model%state_variables)
            write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%state_variables(i)%long_name),trim(model%state_variables(i)%units)
         end do
         do i=1,size(model%bottom_state_variables)
            write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%bottom_state_variables(i)%long_name),trim(model%bottom_state_variables(i)%units)
         end do
         if (add_diagnostic_variables) then
            do i=1,size(model%diagnostic_variables)
               write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%diagnostic_variables(i)%long_name),trim(model%diagnostic_variables(i)%units)
            end do
            do i=1,size(model%horizontal_diagnostic_variables)
               write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%horizontal_diagnostic_variables(i)%long_name),trim(model%horizontal_diagnostic_variables(i)%units)
            end do
         end if
         if (add_conserved_quantities) then
            do i=1,size(model%conserved_quantities)
               write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%conserved_quantities(i)%long_name),trim(model%conserved_quantities(i)%units)
            end do
         end if
         write(out_unit,*)
      case (NETCDF_FMT)
#ifdef NETCDF4
         setn=0
         LEVEL2 'NetCDF version: ',trim(NF90_INQ_LIBVERS())
         iret = nf90_create(output_file,NF90_CLOBBER,ncid)
         call check_err(iret)
!        define dimensions
         iret = nf90_def_dim(ncid, 'lon', 1, lon_dim)
         call check_err(iret)
         iret = nf90_def_dim(ncid, 'lat', 1, lat_dim)
         call check_err(iret)
         iret = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dim)
         call check_err(iret)
         dims(1) = lon_dim; dims(2) = lat_dim; dims(3) = time_dim

!        define coordinates
         iret = nf90_def_var(ncid,'lon',NF90_REAL,lon_dim,lon_id)
         call check_err(iret)
         iret = nf90_def_var(ncid,'lat',NF90_REAL,lat_dim,lat_id)
         call check_err(iret)
         iret = nf90_def_var(ncid,'time',NF90_REAL,time_dim,time_id)
         call check_err(iret)

         write(ncdf_time_str,110) 'seconds',trim(start)
         iret = nf90_put_att(ncid,time_id,"units",trim(ncdf_time_str))
         call check_err(iret)
110 format(A,' since ',A)

!        define variables
         iret = nf90_def_var(ncid,'par',NF90_REAL,dims,par_id)
         call check_err(iret)
         iret = nf90_def_var(ncid,'temp',NF90_REAL,dims,temp_id)
         call check_err(iret)
         iret = nf90_def_var(ncid,'salt',NF90_REAL,dims,salt_id)
         call check_err(iret)

         allocate(statevar_ids(size(model%state_variables)+size(model%bottom_state_variables)))
         allocate(diagnostic_ids(size(model%diagnostic_variables)+size(model%horizontal_diagnostic_variables)))
         allocate(conserved_ids(size(model%conserved_quantities)))
         allocate(conserved_change_ids(size(model%conserved_quantities)))

         do i=1,size(model%state_variables)
            call create_variable(model%state_variables(i),statevar_ids(i))
         end do
         do i=1,size(model%bottom_state_variables)
            call create_variable(model%bottom_state_variables(i),statevar_ids(i+size(model%state_variables)))
         end do

         do i=1,size(model%diagnostic_variables)
            call create_variable(model%diagnostic_variables(i),diagnostic_ids(i))
         end do
         do i=1,size(model%horizontal_diagnostic_variables)
            call create_variable(model%horizontal_diagnostic_variables(i),diagnostic_ids(i+size(model%diagnostic_variables)))
         end do

         do i=1,size(model%conserved_quantities)
            call create_variable(model%conserved_quantities(i),conserved_ids(i))
            iret = nf90_def_var(ncid,'int_'//trim(model%conserved_quantities(i)%name)//'_change',NF90_REAL,dims,conserved_change_ids(i))
            call check_err(iret)
            iret = nf90_put_att(ncid,conserved_change_ids(i),"long_name",'integrated change in '//trim(model%conserved_quantities(i)%long_name))
            call check_err(iret)
            iret = nf90_put_att(ncid,conserved_change_ids(i),"units",trim(model%conserved_quantities(i)%units))
            call check_err(iret)
         end do

!        global attributes
         iret = nf90_put_att(ncid,NF90_GLOBAL,'title',trim(title))
!         history = 'Created by GOTM v. '//RELEASE
!         iret = nf90_put_att(ncid,NF90_GLOBAL,'history',trim(history))
         iret = nf90_put_att(ncid,NF90_GLOBAL,'Conventions','COARDS')
         call check_err(iret)

!        leave define mode
         iret = nf90_enddef(ncid)
         call check_err(iret)

!        save latitude and logitude
         iret = nf90_put_var(ncid,lon_id,longitude)
         iret = nf90_put_var(ncid,lat_id,latitude)

         iret = nf90_sync(ncid)
         call check_err(iret)
#endif
      case default
   end select

   return
96 FATAL 'I could not open ',trim(output_file)
   stop 'init_output'
100 format (A, A, ' (', A, ')')

#ifdef NETCDF4
   contains

   subroutine create_variable(variable,id)
      class (type_external_variable),intent(in)  :: variable
      integer,                       intent(out) :: id
      iret = nf90_def_var(ncid,trim(variable%name),NF90_REAL,dims,id)
      call check_err(iret)
      iret = nf90_put_att(ncid,id,"long_name",trim(variable%long_name))
      call check_err(iret)
      iret = nf90_put_att(ncid,id,"units",trim(variable%units))
      call check_err(iret)
   end subroutine
#endif

   end subroutine init_output
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: do the output
!
! !INTERFACE:
   subroutine do_output(output_format,n)
!
! !DESCRIPTION:
! TODO
   implicit none
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: output_format
   integer(timestepkind), intent(in)   :: n
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL PARAMETERS:
   integer         :: i,j,iret
   integer         :: start(3),edges(3)
   REALTYPE        :: x(1)
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (output_format)
      case(ASCII_FMT)
         call write_time_string(julianday,secondsofday,timestr)
         write (out_unit,FMT='(A)',ADVANCE='NO') timestr
         if (add_environment) then
            write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,par
            write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,temp
            write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,salt
         end if
         do i=1,(size(model%state_variables)+size(model%bottom_state_variables))
            write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,cc(i,1)
         end do
         if (add_diagnostic_variables) then
            do i=1,size(model%diagnostic_variables)
               write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,fabm_get_bulk_diagnostic_data(model,i)
            end do
            do i=1,size(model%horizontal_diagnostic_variables)
               write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,fabm_get_horizontal_diagnostic_data(model,i)
            end do
         end if
         if (add_conserved_quantities) then
            call fabm_get_conserved_quantities(model,totals)
            do i=1,size(model%conserved_quantities)
               write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,totals(i)
            end do
         end if
         write (out_unit,*)

      case (NETCDF_FMT)
#ifdef NETCDF4
         setn = setn + 1
         start(1) = setn; edges(1) = 1

         x(1) = n * timestep
         iret = nf90_put_var(ncid,time_id,x,start,edges)

         start(1) = 1; start(2) = 1; start(3) = setn
         edges(1) = 1; edges(2) = 1; edges(3) = 1

         x(1) = par_sf
         iret = nf90_put_var(ncid,par_id,x,start,edges)
         call check_err(iret)
         x(1) = temp
         iret = nf90_put_var(ncid,temp_id,x,start,edges)
         call check_err(iret)
         x(1) = salt
         iret = nf90_put_var(ncid,salt_id,x,start,edges)
         call check_err(iret)

         do i=1,size(model%state_variables)
            x(1) = cc(i,1)
            iret = nf90_put_var(ncid,statevar_ids(i),x,start,edges)
            call check_err(iret)
         end do

         j=size(model%state_variables)
         do i=1,size(model%bottom_state_variables)
            x(1) = cc(i+j,1)
            iret = nf90_put_var(ncid,statevar_ids(i+j),x,start,edges)
            call check_err(iret)
         end do

         do i=1,size(model%diagnostic_variables)
            x(1) = fabm_get_bulk_diagnostic_data(model,i)
            iret = nf90_put_var(ncid,diagnostic_ids(i),x,start,edges)
            call check_err(iret)
         end do

         j = size(model%diagnostic_variables)
         do i=1,size(model%horizontal_diagnostic_variables)
            x(1) = fabm_get_horizontal_diagnostic_data(model,i)
            iret = nf90_put_var(ncid,diagnostic_ids(i+j),x,start,edges)
            call check_err(iret)
         end do

         call fabm_get_conserved_quantities(model,totals)
         if (n==0_timestepkind) totals0 = totals
         do i=1,size(model%conserved_quantities)
            x(1) = totals(i)
            iret = nf90_put_var(ncid,conserved_ids(i),x,start,edges)
            call check_err(iret)

            x(1) = totals(i)-totals0(i)
            iret = nf90_put_var(ncid,conserved_change_ids(i),x,start,edges)
            call check_err(iret)
         end do
#endif
      case default
   end select

   end subroutine do_output
!EOC

   subroutine check_fabm_expressions()
      class (type_expression),pointer :: expression
      integer :: n
      
      n = 0
      expression => model%root%first_expression
      do while (associated(expression))
         select type (expression)
            class is (type_vertical_integral)
                n = n + 1
         end select
         expression => expression%next
      end do

      allocate(expression_data(n))
      expression_data = _ZERO_

      n = 0
      expression => model%root%first_expression
      do while (associated(expression))
         select type (expression)
            class is (type_vertical_integral)
               n = n + 1
               call fabm_link_horizontal_data(model,trim(expression%output_name),expression_data(n))
         end select
         expression => expression%next
      end do
   end subroutine

   subroutine update_fabm_expressions()
      class (type_expression),pointer :: expression

      expression => model%root%first_expression
      do while (associated(expression))
         select type (expression)
            class is (type_vertical_integral)
               expression%out%p = expression%in%p
               if (.not.expression%average) expression%out%p = expression%out%p*(min(expression%maximum_depth,column_depth)-expression%minimum_depth)
         end select
         expression => expression%next
      end do
   end subroutine
   
#ifdef NETCDF4
   subroutine check_err(iret)
      use netcdf

      integer :: iret

      if (iret .ne. NF90_NOERR) then
         print *, nf90_strerror(iret)
         stop
      endif

   end subroutine
#endif

!-----------------------------------------------------------------------

   end module fabm_0d

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

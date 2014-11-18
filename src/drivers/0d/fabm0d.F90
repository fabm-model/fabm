#include "fabm_driver.h"
#include "fabm_0d.h"
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
   use eqstate,only:rho_feistel

   use fabm
   use fabm_driver
   use fabm_types
   use fabm_expressions
   use fabm_config

   use shared
   use output

   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_run, time_loop, clean_up
!
! !DEFINED PARAMETERS:
   integer, parameter        :: namlst=10, yaml_unit=23
   integer, parameter        :: CENTER=0,SURFACE=1,BOTTOM=2
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!  private data members initialised via namelists
   real(rk)                  :: dt

!  FABM yalmm configuration file
   character(len=PATH_MAX)   :: fabm_yaml_file='fabm.yaml'

   ! Bio model info
   integer  :: ode_method
   logical  :: repair_state
   integer  :: swr_method
   real(rk) :: cloud
   real(rk) :: par_fraction
   real(rk) :: par_background_extinction
   logical  :: apply_self_shading

   ! Environment
   real(rk),target :: current_depth,dens,decimal_yearday
   real(rk)        :: swr_sf,par_sf,par_bt,par_ct,extinction

   real(rk),allocatable :: expression_data(:)

   type (type_bulk_variable_id), save :: id_dens, id_par
   logical                            :: compute_density

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

#define _ODE_ZEROD_
#include "ode_solvers_template.F90"

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Parse the command line
!
! !INTERFACE:
   subroutine  cmdline
   implicit none

!   character(len=*), parameter :: version = '1.0'
   character(len=32) :: arg
   integer :: i
!EOP
!-----------------------------------------------------------------------
!BOC
   do i = 1, command_argument_count()
      call get_command_argument(i, arg)

      select case (arg)
#if 0
      case ('-v', '--version')
         print '(2a)', 'fabm0d version ', version
         stop
#endif
      case ('-h', '--help')
         call print_help()
         stop
      case ('-y', '--yaml')
         call get_command_argument(i+1, fabm_yaml_file)
!         print '(a)', fabm_yaml_file
      case default
#if 0
         print '(a,a,/)', 'Unrecognized command-line option: ', arg
         call print_help()
         stop
#endif
      end select
   end do

   contains

   subroutine print_help()
      print '(a)', 'usage: fabm0d [OPTIONS]'
      print '(a)', ''
      print '(a)', 'Without further options, fabm0d run using default input filenames.'
      print '(a)', ''
      print '(a)', 'fabm0d options:'
      print '(a)', ''
      print '(a)', '  -h, --help        print usage information and exit'
      print '(a)', '  -y, --yaml file   use <file> as FABM yaml configuration file - default fabm.yaml'
   end subroutine print_help

   end subroutine  cmdline

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
   character(len=PATH_MAX)   :: env_file
   integer                   :: i
   real(rk)                  :: depth
   real(rk),parameter        :: invalid_latitude = -100._rk,invalid_longitude = -400.0_rk
   logical                   :: file_exists
   real(rk),allocatable      :: rhs(:,:)

   namelist /model_setup/ title,start,stop,dt,ode_method,repair_state
   namelist /environment/ env_file,swr_method, &
                          latitude,longitude,cloud,par_fraction, &
                          depth,par_background_extinction,apply_self_shading
!EOP
!-----------------------------------------------------------------------
!BOC
   call cmdline

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

   call configure_output(namlst)

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
   column_depth = depth ! Provided depth is the column depth. The modelled biogeochemistry will be positioned at half this depth.
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
   LEVEL1 'init environment'
   LEVEL2 'reading local environment data from:'
   LEVEL3 trim(env_file)
   call init_input()
   call register_input_0d(env_file,1,swr_sf,'shortwave radiation')
   call register_input_0d(env_file,2,temp,'temperature')
   call register_input_0d(env_file,3,salt,'salinity')

   ! Build FABM model tree. Use 'fabm_yaml_file' if available, otherwise fall back to fabm.nml.
   LEVEL1 'initialize FABM'
   LEVEL2 'reading configuration from:'
   inquire(file=trim(fabm_yaml_file),exist=file_exists)
   if (file_exists) then
      ! From YAML file fabm.yaml
      LEVEL3 trim(fabm_yaml_file)
      allocate(model)
      call fabm_create_model_from_yaml_file(model,path=trim(fabm_yaml_file))
   else
      ! From namelists in fabm.nml
      LEVEL3 'fabm.nml'
      model => fabm_create_model_from_file(namlst)
   end if

   ! Send information on spatial domain to FABM (this also allocates memory for diagnostics)
   call fabm_set_domain(model,dt)

   ! Create state variable vector, using the initial values specified by the model,
   ! and link state data to FABM.
   allocate(cc(size(model%state_variables)+size(model%bottom_state_variables)+size(model%surface_state_variables)))
   do i=1,size(model%state_variables)
      cc(i) = model%state_variables(i)%initial_value
      call fabm_link_bulk_state_data(model,i,cc(i))
   end do

   ! Create bottom-bound state variable vector, using the initial values specified by the model,
   ! and link state data to FABM.
   do i=1,size(model%bottom_state_variables)
      cc(size(model%state_variables)+i) = model%bottom_state_variables(i)%initial_value
      call fabm_link_bottom_state_data(model,i,cc(size(model%state_variables)+i))
   end do

   ! Create surface-bound state variable vector, using the initial values specified by the model,
   ! and link state data to FABM.
   do i=1,size(model%surface_state_variables)
      cc(size(model%state_variables)+size(model%bottom_state_variables)+i) = model%surface_state_variables(i)%initial_value
      call fabm_link_surface_state_data(model,i,cc(size(model%state_variables)+size(model%bottom_state_variables)+i))
   end do

   id_dens = fabm_get_bulk_variable_id(model,standard_variables%density)
   compute_density = fabm_variable_needs_values(model,id_dens)
   if (compute_density) call fabm_link_bulk_data(model,id_dens,dens)

   id_par = fabm_get_bulk_variable_id(model,standard_variables%downwelling_photosynthetic_radiative_flux)

   ! Link environmental data to FABM
   call fabm_link_bulk_data(model,standard_variables%temperature,temp)
   call fabm_link_bulk_data(model,standard_variables%practical_salinity,salt)
   if (fabm_variable_needs_values(model,id_par)) call fabm_link_bulk_data(model,id_par,par)
   call fabm_link_bulk_data(model,standard_variables%pressure,current_depth)
   call fabm_link_bulk_data(model,standard_variables%cell_thickness,column_depth)
   call fabm_link_bulk_data(model,standard_variables%depth,current_depth)
   call fabm_link_bulk_data(model,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,extinction)
   call fabm_link_horizontal_data(model,standard_variables%surface_downwelling_photosynthetic_radiative_flux,par_sf)
   call fabm_link_horizontal_data(model,standard_variables%surface_downwelling_shortwave_flux,swr_sf)
   call fabm_link_horizontal_data(model,standard_variables%cloud_area_fraction,cloud)
   call fabm_link_horizontal_data(model,standard_variables%bottom_depth,column_depth)
   call fabm_link_horizontal_data(model,standard_variables%bottom_depth_below_geoid,column_depth)
   if (latitude /=invalid_latitude ) call fabm_link_horizontal_data(model,standard_variables%latitude,latitude)
   if (longitude/=invalid_longitude) call fabm_link_horizontal_data(model,standard_variables%longitude,longitude)
   call fabm_link_scalar_data(model,standard_variables%number_of_days_since_start_of_the_year,decimal_yearday)

   ! Read forcing data specified in input.yaml.
   call init_yaml()

   ! Check whether all dependencies of biogeochemical models have now been fulfilled.
   call fabm_check_ready(model)

   ! Update time and all time-dependent inputs.
   call start_time_step(0_timestepkind)

   ! Allow the model to compute all diagnostics, so output for initial time contains sensible values.
   allocate(rhs(size(cc,1),0:1))
   call get_rhs(.true.,size(cc,1),cc,rhs)

   ! Output variable values at initial time
   LEVEL1 'init_output'
   call init_output(start)
   call do_output(0_timestepkind)

   STDERR LINE

   return

90 FATAL 'I could not open run.nml for reading'
   stop 'init_run'
91 FATAL 'I could not read the "model_setup" namelist'
   stop 'init_run'
92 FATAL 'I could not read the "environment" namelist'
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
      if (associated(pair)) call driver%log_message('Forcing data specified in input.yaml:')
      do while (associated(pair))
         variable_name = trim(pair%key)
         if (variable_name=='') call driver%fatal_error('init_yaml_input','Empty variable name specified.')
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

      type (type_error),         pointer :: config_error
      class (type_node),         pointer :: node
      class (type_scalar),       pointer :: scalar
      character(len=1024)                :: path,message
      integer                            :: column
      real(rk)                           :: scale_factor
      real(rk)                           :: relaxation_time
      logical                            :: is_state_variable
      type (type_input_data),    pointer :: input_data
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

      ! Prepend to list of input data.
      allocate(input_data)
      input_data%next => first_input_data
      first_input_data => input_data
      input_data%variable_name = variable_name
      call driver%log_message('  '//trim(input_data%variable_name))

      scalar => mapping%get_scalar('constant_value',required=.false.,error=config_error)
      if (associated(scalar)) then
         ! Input variable is set to a constant value.
         input_data%value = mapping%get_real('constant_value',error=config_error)
         if (associated(config_error)) call fatal_error('parse_input_variable',config_error%message)

         ! Make sure keys related to time-varying input are not present.
         node => mapping%get('file')
         if (associated(node)) call fatal_error('parse_input_variable','input.yaml, variable "'//trim(variable_name)//'": keys "constant_value" and "file" cannot both be present.')
         node => mapping%get('column')
         if (associated(node)) call fatal_error('parse_input_variable','input.yaml, variable "'//trim(variable_name)//'": keys "constant_value" and "column" cannot both be present.')
         node => mapping%get('scale_factor')
         if (associated(node)) call fatal_error('parse_input_variable','input.yaml, variable "'//trim(variable_name)//'": keys "constant_value" and "scale_factor" cannot both be present.')
         write (message,'(g13.6)') input_data%value
         call driver%log_message('    constant_value = '//adjustl(message))
      else
         ! Input variable is set to a time-varying value. Obtain path, column number and scale factor.
         path = mapping%get_string('file',error=config_error)
         if (associated(config_error)) call fatal_error('parse_input_variable',config_error%message)
         column = mapping%get_integer('column',default=1,error=config_error)
         if (associated(config_error)) call fatal_error('parse_input_variable',config_error%message)
         scale_factor = mapping%get_real('scale_factor',default=1.0_rk,error=config_error)
         if (associated(config_error)) call fatal_error('parse_input_variable',config_error%message)
         call register_input_0d(path,column,input_data%value,variable_name,scale_factor=scale_factor)
         call driver%log_message('    file = '//trim(path))
         write (message,'(i0)') column
         call driver%log_message('    column = '//adjustl(message))
         write (message,'(g13.6)') scale_factor
         call driver%log_message('    scale factor = '//adjustl(message))
      end if

      if (is_state_variable) then
         ! This is a state variable. Obtain associated relaxation time.
         relaxation_time = mapping%get_real('relaxation_time',default=1.e15_rk,error=config_error)
         if (associated(config_error)) call fatal_error('parse_input_variable',config_error%message)
      else
         ! This is not a state variable. Make sure no relaxation time is specified.
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

   subroutine start_time_step(n)
      integer(timestepkind),intent(in) :: n

      real(rk)                         :: bio_albedo

      ! Update time in time manager
      call update_time(n)

      ! Compute decimal year day (input for some biogeochemical models)
      decimal_yearday = yearday-1 + dble(secondsofday)/86400

      ! Update environment (i.e., read from input files)
      call do_input(julianday,secondsofday)

      ! Calculate light extinction
      extinction = 0.0_rk
      if (apply_self_shading) call fabm_get_light_extinction(model,extinction)
      extinction = extinction + par_background_extinction

      ! Calculate photosynthetically active radiation at surface, if it is not provided in the input file.
      if (swr_method==0) then
         ! Calculate photosynthetically active radiation from geographic location, time, cloud cover.
         call fabm_get_albedo(model,bio_albedo)
         swr_sf = short_wave_radiation(julianday,secondsofday,longitude,latitude,cloud,bio_albedo)
      end if

      ! Multiply by fraction of short-wave radiation that is photosynthetically active.
      par_sf = par_fraction*swr_sf

      ! Apply light attentuation with depth, unless local light is provided in the input file.
      if (swr_method/=2) then
         ! Either we calculate surface PAR, or surface PAR is provided.
         ! Calculate the local PAR at the given depth from par fraction, extinction coefficient, and depth.
         par_ct = par_sf*exp(-0.5_rk*column_depth*extinction)
         par_bt = par_sf*exp(-column_depth*extinction)
      else
         par_ct = par_sf
         par_bt = par_sf
      end if
      call update_depth(CENTER)

      call fabm_get_light(model)

      ! Compute density from temperature and salinity, if required by biogeochemistry.
      if (compute_density) dens = rho_feistel(salt,temp,5._rk*column_depth,.true.)
   end subroutine start_time_step

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
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL1 'time_loop'

   do n=MinN,MaxN
      ! Update time and all time-dependent inputs.
      call start_time_step(n)

      ! Repair state before calling FABM
      call do_repair_state('0d::time_loop(), before ode_solver()')

      call fabm_update_time(model,real(n,rk))

      ! Integrate one time step
      call ode_solver(ode_method,size(model%state_variables)+size(model%bottom_state_variables),dt,cc,get_rhs,get_ppdd)

      ! ODE solver may have redirected the current state with to an array with intermediate values.
      ! Reset to global array.
      do i=1,size(model%state_variables)
         call fabm_link_bulk_state_data(model,i,cc(i))
      end do
      do i=1,size(model%bottom_state_variables)
         call fabm_link_bottom_state_data(model,i,cc(size(model%state_variables)+i))
      end do
      do i=1,size(model%surface_state_variables)
         call fabm_link_surface_state_data(model,i,cc(size(model%state_variables)+size(model%bottom_state_variables)+i))
      end do

      call do_repair_state('0d::time_loop(), after ode_solver()')

      ! Do output
      call do_output(n)

   end do
   STDERR LINE

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
   call clean_output()

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
      stop 'fabm0d::do_repair_state'
   end if

   end subroutine do_repair_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the right-hand side of the ODE system.
!
! !INTERFACE:
   subroutine get_ppdd(first,numc,cc,pp,dd)
!
! !DESCRIPTION:
! TODO
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: numc
   real(rk), intent(in)                 :: cc(1:numc)
!
! !OUTPUT PARAMETERS:
   real(rk), intent(out)                :: pp(1:numc,1:numc)
   real(rk), intent(out)                :: dd(1:numc,1:numc)
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
   pp = 0.0_rk
   dd = 0.0_rk

   ! Send current state to FABM
   ! (this may differ from the global state cc if using a multi-step integration scheme such as Runge-Kutta)
   do n=1,size(model%state_variables)
      call fabm_link_bulk_state_data(model,n,cc(n))
   end do
   do n=1,size(model%bottom_state_variables)
      call fabm_link_bottom_state_data(model,n,cc(size(model%state_variables)+n))
   end do
   do n=1,size(model%surface_state_variables)
      call fabm_link_surface_state_data(model,n,cc(size(model%state_variables)+size(model%bottom_state_variables)+n))
   end do

   ! Shortcut to the number of pelagic state variables.
   n = size(model%state_variables)

   ! Calculate temporal derivatives due to benthic processes.
   call update_depth(BOTTOM)
   call fabm_do_benthos(model,pp,dd,n)

   ! For pelagic variables: translate bottom flux to into change in concentration
   pp(1:n,:) = pp(1:n,:)/column_depth
   dd(1:n,:) = dd(1:n,:)/column_depth

   ! For pelagic variables: surface and bottom flux (rate per surface area) to concentration (rate per volume)
   call update_depth(CENTER)
   call fabm_do(model,pp,dd)

   end subroutine get_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the right-hand side of the ODE system.
!
! !INTERFACE:
   subroutine get_rhs(first,numc,cc,rhs)
!
! !DESCRIPTION:
! TODO
!
! !INPUT PARAMETERS:
   logical, intent(in)                  :: first
   integer, intent(in)                  :: numc
   real(rk), intent(in)                 :: cc(1:numc)
!
! !OUTPUT PARAMETERS:
   real(rk), intent(out)                :: rhs(1:numc)
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
   rhs = 0.0_rk

   ! Send current state to FABM
   ! (this may differ from the global state cc if using a multi-step integration scheme such as Runge-Kutta)
   do n=1,size(model%state_variables)
      call fabm_link_bulk_state_data(model,n,cc(n))
   end do
   do n=1,size(model%bottom_state_variables)
      call fabm_link_bottom_state_data(model,n,cc(size(model%state_variables)+n))
   end do
   do n=1,size(model%surface_state_variables)
      call fabm_link_surface_state_data(model,n,cc(size(model%state_variables)+size(model%bottom_state_variables)+n))
   end do

   ! Shortcut to the number of pelagic state variables.
   n = size(model%state_variables)

   ! Calculate temporal derivatives due to surface-bound processes.
   call update_depth(SURFACE)
   call fabm_do_surface(model,rhs(1:n),rhs(n+size(model%bottom_state_variables)+1:))

   ! Calculate temporal derivatives due to bottom-bound processes.
   call update_depth(BOTTOM)
   call fabm_do_bottom(model,rhs(1:n),rhs(n+1:n+size(model%bottom_state_variables)))

   ! For pelagic variables: surface and bottom flux (rate per surface area) to concentration (rate per volume)
   rhs(1:n) = rhs(1:n)/column_depth

   ! Add change in pelagic variables.
   call update_depth(CENTER)
   call fabm_do(model,rhs)

   end subroutine get_rhs
!EOC

!-----------------------------------------------------------------------

   end module fabm_0d

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#include "fabm_driver.h"
#include "fabm_0d.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Output manager
!
! !INTERFACE:
   module output
!
! !DESCRIPTION:
! TODO
!
   ! From FABM
   use fabm
   use fabm_types
   use fabm_driver
   use fabm_properties

   ! From GOTM
   use time
#ifdef NETCDF4
   use field_manager
   use output_manager_core, only:output_manager_host=>host, type_output_manager_host=>type_host
   use output_manager
#endif

   use shared

   implicit none

   private

   public configure_output, init_output, do_output, clean_output

#ifdef NETCDF4
   public register_output_fields, fm

   type,extends(type_output_manager_host) :: type_fabm0d_host
   contains
      procedure :: julian_day    => fabm0d_host_julian_day
      procedure :: calendar_date => fabm0d_host_calendar_date
   end type

   type (type_field_manager), target :: fm
#endif

   integer, parameter :: ASCII_FMT  = 1
   integer, parameter :: NETCDF_FMT = 2

   character(len=PATH_MAX) :: output_file
   integer, public :: output_format
   logical :: add_environment
   logical :: add_conserved_quantities
   logical :: add_diagnostic_variables
   integer(timestepkind) :: nsave

   integer                    :: out_unit = -1

   character, parameter       :: separator = char(9)
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: configure output from namelists or YAML
!
! !INTERFACE:
   subroutine configure_output(namlst)
!
! !DESCRIPTION:
! TODO
!
! !INPUT PARAMETERS:
   integer, intent(in) :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL PARAMETERS:
   namelist /output/ output_file,output_format,nsave,add_environment, &
                     add_diagnostic_variables, add_conserved_quantities
   integer :: ios
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read output namelist
   output_file = ''
   output_format = ASCII_FMT
   nsave = 1
   add_environment = .false.
   add_conserved_quantities = .false.
   add_diagnostic_variables = .false.

   read(namlst,nml=output,iostat=ios)
   if (ios/=0) call driver%fatal_error('configure_output','run.nml: I could not read the "output" namelist.')
   if (output_file=='') call driver%fatal_error('configure_output','run.nml: "output_file" must be set to a valid file path in "output" namelist.')

   end subroutine configure_output
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: prepare for output
!
! !INTERFACE:
   subroutine init_output(start)
!
! !DESCRIPTION:
! TODO
!
! !INPUT PARAMETERS:
   character(len=*), intent(in) :: start
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL PARAMETERS:
   integer                        :: i,iret
   type (type_input_data),pointer :: input_data,input_data2
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef NETCDF4  
   if (output_format /= 1) then
      call send_output_data()
      return
   end if
#endif

   LEVEL2 'writing results to:'
   LEVEL3 trim(output_file)
   select case (output_format)
      case (ASCII_FMT)
         out_unit = get_free_unit()
         open(out_unit,file=trim(output_file),action='write',status='replace',err=96)
         LEVEL3 'ASCII format'

         ! Write header to the output file.
         write(out_unit,FMT='(''# '',A)') title
         write(out_unit,FMT='(''# '',A)',ADVANCE='NO') 'time'
         if (add_environment) then
            write(out_unit,FMT=100,ADVANCE='NO') separator,'photosynthetically active radiation','W/m2'
            write(out_unit,FMT=100,ADVANCE='NO') separator,'temperature',                        'degrees C'
            write(out_unit,FMT=100,ADVANCE='NO') separator,'salinity',                           'kg/m3'
         end if
         do i=1,size(model%state_variables)
            if (model%state_variables(i)%output/=output_none) &
               write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%state_variables(i)%long_name), &
                  trim(model%state_variables(i)%units)
         end do
         do i=1,size(model%bottom_state_variables)
            if (model%bottom_state_variables(i)%output/=output_none) &
               write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%bottom_state_variables(i)%long_name), &
                  trim(model%bottom_state_variables(i)%units)
         end do
         do i=1,size(model%surface_state_variables)
            if (model%surface_state_variables(i)%output/=output_none) &
               write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%surface_state_variables(i)%long_name), &
                  trim(model%surface_state_variables(i)%units)
         end do
         if (add_diagnostic_variables) then
            do i=1,size(model%diagnostic_variables)
               if (model%diagnostic_variables(i)%output/=output_none) &
                  write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%diagnostic_variables(i)%long_name), &
                     trim(model%diagnostic_variables(i)%units)
            end do
            do i=1,size(model%horizontal_diagnostic_variables)
               if (model%horizontal_diagnostic_variables(i)%output/=output_none) &
                  write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%horizontal_diagnostic_variables(i)%long_name), &
                     trim(model%horizontal_diagnostic_variables(i)%units)
            end do
         end if
         if (add_conserved_quantities) then
            compute_conserved_quantities = .true.
            do i=1,size(model%conserved_quantities)
               if (model%conserved_quantities(i)%output/=output_none) &
                  write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%conserved_quantities(i)%long_name), &
                     trim(model%conserved_quantities(i)%units)
            end do
         end if
         write(out_unit,*)
      case (NETCDF_FMT)
      case default
   end select

   return
96 call driver%fatal_error('init_output','I could not open '//trim(output_file))

100 format (A, A, ' (', A, ')')

   end subroutine init_output
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: do the output
!
! !INTERFACE:
   subroutine do_output(n)
!
! !DESCRIPTION:
! TODO
!
! !INPUT PARAMETERS:
   integer(timestepkind), intent(in) :: n
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL PARAMETERS:
   integer         :: i,j,iret
   integer         :: start(3)
   type (type_input_data),pointer :: input_data
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef NETCDF4  
   if (output_format /= 1) then
      call output_manager_save(julianday,secondsofday,int(n))
      return
   end if
#endif

   if (mod(n,nsave)/=0) return

   select case (output_format)
      case(ASCII_FMT)
         call write_time_string(julianday,secondsofday,timestr)
         write (out_unit,FMT='(A)',ADVANCE='NO') timestr
         if (add_environment) then
            write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,par
            write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,temp
            write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,salt
         end if
         do i=1,size(model%state_variables)
            if (model%state_variables(i)%output/=output_none) &
               write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,model%get_data(model%state_variables(i)%globalid)
         end do
         do i=1,size(model%bottom_state_variables)
            if (model%bottom_state_variables(i)%output/=output_none) &
               write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,model%get_data(model%bottom_state_variables(i)%globalid)
         end do
         do i=1,size(model%surface_state_variables)
            if (model%surface_state_variables(i)%output/=output_none) &
               write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,model%get_data(model%surface_state_variables(i)%globalid)
         end do
         if (add_diagnostic_variables) then
            do i=1,size(model%diagnostic_variables)
               if (model%diagnostic_variables(i)%output/=output_none) &
                  write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,fabm_get_bulk_diagnostic_data(model,i)
            end do
            do i=1,size(model%horizontal_diagnostic_variables)
               if (model%horizontal_diagnostic_variables(i)%output/=output_none) &
                  write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,fabm_get_horizontal_diagnostic_data(model,i)
            end do
         end if
         if (add_conserved_quantities) then
            do i=1,size(model%conserved_quantities)
               if (model%conserved_quantities(i)%output/=output_none) &
                  write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,totals(i)
            end do
         end if
         write (out_unit,*)

      case (NETCDF_FMT)
      case default
   end select

   end subroutine do_output
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Clean up.
!
! !INTERFACE:
   subroutine clean_output(ignore_errors)
!
! !DESCRIPTION:
! Close all open files.
!
! !INPUT PARAMETERS:
   logical, intent(in) :: ignore_errors
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL PARAMETERS:
   integer :: iret
!EOP
!-----------------------------------------------------------------------
!BOC

   if (out_unit/=-1) close(out_unit)

#ifdef NETCDF4
   call output_manager_clean()
   call fm%finalize()
#endif

   end subroutine clean_output
!EOC

#ifdef NETCDF4
   subroutine register_output_fields()

   integer :: i
   logical :: in_output

   LEVEL1 'output_manager'
   allocate(type_fabm0d_host::output_manager_host)
   call output_manager_init(fm,title)

   LEVEL1 'field_manager'
   call fm%register_dimension('lon',1,id=id_dim_lon)
   call fm%register_dimension('lat',1,id=id_dim_lat)
   call fm%register_dimension('time',id=id_dim_time)
   call fm%initialize(prepend_by_default=(/id_dim_lon,id_dim_lat/),append_by_default=(/id_dim_time/))
   call fm%register('lon','degrees_east','longitude',dimensions=(/id_dim_lon/),no_default_dimensions=.true.,data0d=longitude,coordinate_dimension=id_dim_lon)
   call fm%register('lat','degrees_north','latitude',dimensions=(/id_dim_lat/),no_default_dimensions=.true.,data0d=latitude,coordinate_dimension=id_dim_lat)

   call fm%register('par','W/m^2','par',standard_name='downwelling_photosynthetic_radiative_flux',data0d=par)
   call fm%register('temp','Celsius','temperature',standard_name='sea_water_temperature',data0d=temp)
   call fm%register('salt','1e-3','salinity',standard_name='sea_water_practical_salinity',data0d=salt)

   ! state variables
   do i=1,size(model%state_variables)
      in_output = register(model%state_variables(i))
   end do
   do i=1,size(model%bottom_state_variables)
      in_output = register(model%bottom_state_variables(i))
   end do
   do i=1,size(model%surface_state_variables)
      in_output = register(model%surface_state_variables(i))
   end do

   ! diagnostic variables
   do i=1,size(model%diagnostic_variables)
      model%diagnostic_variables(i)%save = register(model%diagnostic_variables(i))
   end do
   do i=1,size(model%horizontal_diagnostic_variables)
      model%horizontal_diagnostic_variables(i)%save = register(model%horizontal_diagnostic_variables(i))
   end do

   ! conserved quantities
   do i=1,size(model%conserved_quantities)
      call fm%register('int_change_in_'//trim(model%conserved_quantities(i)%name), trim(model%conserved_quantities(i)%units)//'*m', 'integrated change in '//trim(model%conserved_quantities(i)%long_name), &
                        minimum=model%conserved_quantities(i)%minimum, maximum=model%conserved_quantities(i)%maximum, fill_value=model%conserved_quantities(i)%missing_value, &
                        category='fabm/conservation', output_level=output_level_default, used=in_output, data0d=int_change_in_totals(i))
      if (in_output) compute_conserved_quantities = .true.
   end do

   contains

   function register(variable) result(used)
      class (type_external_variable),intent(in) :: variable
      logical                                   :: used

      integer                        :: output_level
      type (type_field),     pointer :: field
      class (type_property), pointer :: property

      output_level = output_level_default
      if (variable%output==output_none) output_level = output_level_debug
      call fm%register(variable%name, variable%units, variable%long_name, &
                       minimum=variable%minimum, maximum=variable%maximum, fill_value=variable%missing_value, &
                       category='fabm'//variable%target%owner%get_path(), output_level=output_level, used=used, field=field)
      property => variable%properties%first
      do while (associated(property))
         select type (property)
         class is (type_real_property)
            call field%set_attribute(property%name,property%value)
         end select
         property => property%next
      end do
   end function register

   end subroutine register_output_fields

   subroutine send_output_data()
      integer          :: i
      real(rk),pointer :: pdata

      do i=1,size(model%state_variables)
         pdata => model%get_data(model%state_variables(i)%globalid)
         call fm%send_data(model%state_variables(i)%name, pdata)
      end do
      do i=1,size(model%bottom_state_variables)
         pdata => model%get_data(model%bottom_state_variables(i)%globalid)
         call fm%send_data(model%bottom_state_variables(i)%name, pdata)
      end do
      do i=1,size(model%surface_state_variables)
         pdata => model%get_data(model%surface_state_variables(i)%globalid)
         call fm%send_data(model%surface_state_variables(i)%name, pdata)
      end do

      do i=1,size(model%diagnostic_variables)
         if (model%diagnostic_variables(i)%save) &
            call fm%send_data(model%diagnostic_variables(i)%name, fabm_get_interior_diagnostic_data(model,i))
      end do
      do i=1,size(model%horizontal_diagnostic_variables)
         if (model%horizontal_diagnostic_variables(i)%save) &
            call fm%send_data(model%horizontal_diagnostic_variables(i)%name, fabm_get_horizontal_diagnostic_data(model,i))
      end do
   end subroutine send_output_data

   subroutine fabm0d_host_julian_day(self,yyyy,mm,dd,julian)
      class (type_fabm0d_host), intent(in) :: self
      integer, intent(in)  :: yyyy,mm,dd
      integer, intent(out) :: julian
      call julian_day(yyyy,mm,dd,julian)
   end subroutine

   subroutine fabm0d_host_calendar_date(self,julian,yyyy,mm,dd)
      class (type_fabm0d_host), intent(in) :: self
      integer, intent(in)  :: julian
      integer, intent(out) :: yyyy,mm,dd
      call calendar_date(julian,yyyy,mm,dd)
   end subroutine
#endif

!-----------------------------------------------------------------------

   end module output

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

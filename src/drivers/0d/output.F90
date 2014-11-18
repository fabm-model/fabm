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
#ifdef NETCDF4
   use netcdf
#endif

   ! From FABM
   use fabm
   use fabm_types

   ! From GOTM
   use time

   use shared

   implicit none

   private

   public configure_output,init_output,do_output,clean_output

   integer, parameter :: ASCII_FMT  = 1
   integer, parameter :: NETCDF_FMT = 2

   character(len=PATH_MAX) :: output_file
   integer :: output_format
   logical :: add_environment
   logical :: add_conserved_quantities
   logical :: add_diagnostic_variables
   integer(timestepkind) :: nsave

   integer                    :: out_unit = -1

   character, parameter       :: separator = char(9)

   real(rk), allocatable, dimension(:) :: totals,totals0,totals_hz

#ifdef NETCDF4
   integer                   :: ncid = -1
   integer                   :: setn
   integer                   :: time_id
   integer                   :: par_id,temp_id,salt_id
   integer, allocatable, dimension(:) :: statevar_ids,diagnostic_ids,conserved_ids,conserved_change_ids
#endif
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

   read(namlst,nml=output,err=93)

   if (output_file=='') then
      FATAL 'run.nml: "output_file" must be set to a valid file path in "output" namelist.'
      stop 'configure_output'
   end if

   return

93 FATAL 'run.nml: I could not read the "output" namelist.'
   stop 'configure_output'

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
#ifdef NETCDF4
   integer         :: lon_dim,lat_dim,time_dim
   integer         :: lon_id,lat_id
   integer         :: dims(3)
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
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
            write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%state_variables(i)%long_name),trim(model%state_variables(i)%units)
         end do
         do i=1,size(model%bottom_state_variables)
            write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%bottom_state_variables(i)%long_name),trim(model%bottom_state_variables(i)%units)
         end do
         do i=1,size(model%surface_state_variables)
            write(out_unit,FMT=100,ADVANCE='NO') separator,trim(model%surface_state_variables(i)%long_name),trim(model%surface_state_variables(i)%units)
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
         LEVEL3 'NetCDF version: ',trim(NF90_INQ_LIBVERS())
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

         iret = nf90_put_att(ncid,time_id,'units','seconds since '//trim(start))
         call check_err(iret)

!        define variables
         input_data => first_input_data
         do while (associated(input_data))
            ! First check if an input variable of this name has already been added to NetCDF.
            ! If so, this predefined one takes precedence, and we do nothing.
            input_data2 => first_input_data
            do while (.not.associated(input_data2,input_data))
               if (input_data2%variable_name==input_data%variable_name) exit
               input_data2 => input_data2%next
            end do

            if (associated(input_data2,input_data)) then
               iret = nf90_def_var(ncid,input_data%variable_name,NF90_REAL,dims,input_data%ncid)
               call check_err(iret)
            end if
            input_data => input_data%next
         end do
         iret = nf90_def_var(ncid,'par',NF90_REAL,dims,par_id)
         call check_err(iret)
         iret = nf90_def_var(ncid,'temp',NF90_REAL,dims,temp_id)
         call check_err(iret)
         iret = nf90_def_var(ncid,'salt',NF90_REAL,dims,salt_id)
         call check_err(iret)

         allocate(statevar_ids(size(model%state_variables)+size(model%bottom_state_variables)+size(model%surface_state_variables)))
         allocate(diagnostic_ids(size(model%diagnostic_variables)+size(model%horizontal_diagnostic_variables)))
         allocate(conserved_ids(size(model%conserved_quantities)))
         allocate(conserved_change_ids(size(model%conserved_quantities)))

         do i=1,size(model%state_variables)
            call create_variable(model%state_variables(i),statevar_ids(i))
         end do
         do i=1,size(model%bottom_state_variables)
            call create_variable(model%bottom_state_variables(i),statevar_ids(i+size(model%state_variables)))
         end do
         do i=1,size(model%surface_state_variables)
            call create_variable(model%surface_state_variables(i),statevar_ids(i+size(model%state_variables)+size(model%bottom_state_variables)))
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
            iret = nf90_put_att(ncid,conserved_change_ids(i),'long_name','integrated change in '//trim(model%conserved_quantities(i)%long_name))
            call check_err(iret)
            iret = nf90_put_att(ncid,conserved_change_ids(i),'units',trim(model%conserved_quantities(i)%units))
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

   ! Allocate space for totals of conserved quantities.
   allocate(totals0  (size(model%conserved_quantities)))  ! at initial time
   allocate(totals   (size(model%conserved_quantities)))  ! at current time
   allocate(totals_hz(size(model%conserved_quantities)))  ! at current time, on top/bottom interfaces only

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
      iret = nf90_put_att(ncid,id,'long_name',trim(variable%long_name))
      call check_err(iret)
      iret = nf90_put_att(ncid,id,'units',trim(variable%units))
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
         do i=1,size(cc)
            write (out_unit,FMT='(A,E16.8E3)',ADVANCE='NO') separator,cc(i)
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
         start(1) = setn

         iret = nf90_put_var(ncid,time_id,n*timestep,start(1:1))

         start(1) = 1; start(2) = 1; start(3) = setn

         input_data => first_input_data
         do while (associated(input_data))
            if (input_data%ncid/=-1) then
               iret = nf90_put_var(ncid,input_data%ncid,input_data%value,start)
               call check_err(iret)
            end if
            input_data => input_data%next
         end do
         iret = nf90_put_var(ncid,par_id,par,start)
         call check_err(iret)
         iret = nf90_put_var(ncid,temp_id,temp,start)
         call check_err(iret)
         iret = nf90_put_var(ncid,salt_id,salt,start)
         call check_err(iret)

         do i=1,size(cc)
            iret = nf90_put_var(ncid,statevar_ids(i),cc(i),start)
            call check_err(iret)
         end do

         do i=1,size(model%diagnostic_variables)
            iret = nf90_put_var(ncid,diagnostic_ids(i),fabm_get_bulk_diagnostic_data(model,i),start)
            call check_err(iret)
         end do

         j = size(model%diagnostic_variables)
         do i=1,size(model%horizontal_diagnostic_variables)
            iret = nf90_put_var(ncid,diagnostic_ids(i+j),fabm_get_horizontal_diagnostic_data(model,i),start)
            call check_err(iret)
         end do

         call fabm_get_conserved_quantities(model,totals)
         call fabm_get_horizontal_conserved_quantities(model,totals_hz)
         totals = totals + totals_hz/column_depth
         if (n==0_timestepkind) totals0 = totals
         do i=1,size(model%conserved_quantities)
            iret = nf90_put_var(ncid,conserved_ids(i),totals(i),start)
            call check_err(iret)

            iret = nf90_put_var(ncid,conserved_change_ids(i),totals(i)-totals0(i),start)
            call check_err(iret)
         end do
#endif
      case default
   end select

   end subroutine do_output
!EOC

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
!BOP
!
! !IROUTINE: Clean up.
!
! !INTERFACE:
   subroutine clean_output()
!
! !DESCRIPTION:
! Close all open files.
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
   if (ncid/=-1) then
      iret = nf90_close(ncid)
      call check_err(iret)
   end if
#endif

   end subroutine clean_output
!EOC

!-----------------------------------------------------------------------

   end module output

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

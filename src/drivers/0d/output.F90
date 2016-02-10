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

   use shared

   implicit none

   private

   public configure_output,init_output,do_output,clean_output

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

   real(rk), allocatable, dimension(:) :: totals,totals0,totals_hz

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

   if (output_file=='') call driver%fatal_error('configure_output','run.nml: "output_file" must be set to a valid file path in "output" namelist.')

   return

93 call driver%fatal_error('configure_output','run.nml: I could not read the "output" namelist.')

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

   ! Allocate space for totals of conserved quantities.
   allocate(totals0  (size(model%conserved_quantities)))  ! at initial time
   allocate(totals   (size(model%conserved_quantities)))  ! at current time
   allocate(totals_hz(size(model%conserved_quantities)))  ! at current time, on top/bottom interfaces only

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
            call fabm_get_conserved_quantities(model,totals)
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

   end subroutine clean_output
!EOC

!-----------------------------------------------------------------------

   end module output

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

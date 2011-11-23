!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ocean_tpm_mod  !{
!
!<CONTACT EMAIL="Richard.Slater@noaa.gov"> Richard D. Slater
!</CONTACT>
!
!<REVIEWER EMAIL="John.Dunne@noaa.gov"> John P. Dunne
!</REVIEWER>
!
!<OVERVIEW>
! Ocean tracer package module
!</OVERVIEW>
!
!<DESCRIPTION>
!       Currently this module only works for the ocean model,
!       but it could be extended (or generalized) to work with other
!       models.
!
!       This module consists of eight subroutines, three are called as
!       the model is intialized, four are called every time-step, and
!       one is called at model ending. The subroutines are called in
!       the following order.
!
!       These routines are called once at model startup in the
!       ocean_tracer_init routine:
!
!       ocean_tpm_init: This routine saves pointers to "global" model
!               structures, such as Grid and Domain. Also this
!               routine will call specified routines to set default
!               values for each tracer for such things as advection
!               scheme, tracer name, etc.
!
!       ocean_tpm_start: This routine calls specified routines to
!               allocate appropriate storage for the tracer packages,
!               perform pre-processing and initialization (possibly
!               from extra restart information) and set parameters,
!               either via namelist or via the field manager.
!
!       These routines are called each time-step from
!       update_ocean_tracer (one before integration and one after):
!
!       ocean_tpm_sbc: Calls specified routines to handle surface
!               coundary condition calculations. Some or all of
!               this functionality may be moved into a new, generalized
!               boundary condition manager.
!
!       ocean_tpm_bbc: Calls specified routines to handle bottom
!               coundary condition calculations.
!
!       ocean_tpm_source: Calls specified routines to calculate the
!               source array for each tracer in the tracer packages.
!
!       ocean_tpm_tracer: For those packages which need to do
!               post-processing after the continuity equation has
!               been integrated, calls may be placed here. This
!               could be for global, annual means, for instance.
!
!       This routine is called once at the end of the run from
!       ocean_tracer_end:
!
!       ocean_tpm_end: Call routines to finish up any loose ends, such
!               as saving extra restart fields.
!</DESCRIPTION>
!
! <INFO>
! </INFO>
!
! $Id: ocean_tpm.F90 66 2010-03-12 20:16:19Z jornbr $
!

use field_manager_mod
use mpp_mod,            only: stdout, mpp_error, FATAL
use ocean_tpm_util_mod, only: dtts
use ocean_types_mod,    only: ocean_thickness_type, ocean_density_type

!
!       Place tracer modules here
!

use ocean_age_tracer_mod
use ocmip2_cfc_mod
use ocmip2_biotic_mod
use ocean_bgc_restore_mod
use ocean_bgc_phyto_mod
use ocean_tracer_fabm_mod  ! JornB

!
!       force all variables to be "typed"
!

implicit none

!
!       Set all variables to be private by default

private

!
!       Private routines
!

private do_time_calc

!
!       Public routines
!

public ocean_tpm_bbc
public ocean_tpm_end
public ocean_tpm_init
public ocean_tpm_sbc
public ocean_tpm_source
public ocean_tpm_start
public ocean_tpm_tracer

!
!       private parameters
!

character(len=48), parameter                    :: mod_name = 'ocean_tpm_mod'

!
!       Public variables
!
!
!       Private variables
!

contains

!#######################################################################
! <SUBROUTINE NAME="do_time_calc">
!
! <DESCRIPTION>
!       call subroutines to perform time calculations
! </DESCRIPTION>
!

subroutine do_time_calc  !{

!
!-----------------------------------------------------------------------
!
!       Modules
!
!-----------------------------------------------------------------------
!

use time_manager_mod, only: time_type, set_time, get_time,      &
     set_date, get_date, days_in_month, operator(+),            &
     operator(.le.), operator(.eq.),                            &
     operator(*), operator(/), operator(-)
use ocean_tpm_util_mod, only: time, taum1, tau, taup1,          &
     imonth, iyear, end_of_day, end_of_month, end_of_year, mid_month

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'do_time_calc'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       Local variables
!-----------------------------------------------------------------------
!

integer         :: length
type(time_type) :: target_time
type(time_type) :: temp_time
type(time_type) :: dt_time
integer         :: isec
integer         :: iday
real            :: dayint
integer         :: days
integer         :: months
integer         :: years
integer         :: hours
integer         :: minutes
integer         :: seconds
integer, save   :: time_tau = -1000
integer         :: isec2
integer         :: iday2
real            :: daymodel

!
!-----------------------------------------------------------------------
!       Return if this routine has already been called this time-step
!-----------------------------------------------------------------------
!

if (time_tau .eq. time%tau) then  !{
  return
endif  !}

!
!-----------------------------------------------------------------------
!       Set up some things
!-----------------------------------------------------------------------
!

taum1 = time%taum1
tau   = time%tau
taup1 = time%taup1

time_tau = tau

dt_time = set_time (seconds=int(dtts), days=0)

call get_date (time%model_time, years, months, days,            &
               hours, minutes, seconds)

!
!-----------------------------------------------------------------------
!     is it within 1/2 time step of the end of a day ?
!-----------------------------------------------------------------------
!

end_of_day = set_switch (1.0, time%model_time, dt_time)

!
!-----------------------------------------------------------------------
!     is it within 1/2 time step of the middle of a month ?
!-----------------------------------------------------------------------
!

length = days_in_month(time%model_time)
temp_time = set_time(0, length)/2
target_time = set_date(years, months, 1) + temp_time
call get_time (target_time, isec, iday)
dayint = iday + isec/86400.0
mid_month = set_switch (dayint, time%model_time, dt_time)

!
!-----------------------------------------------------------------------
!     is it within 1/2 time step of the end of the month ?
!-----------------------------------------------------------------------
!

length = days_in_month(time%model_time)
target_time = set_date(years, months, 1) + set_time(0, length)
call get_time (target_time, isec, iday)
dayint = iday + isec/86400.
if (days .eq. 1 .and. hours .eq. 0 .and. minutes .eq. 0 .and.   &
    seconds .eq. 0) dayint = dayint - length
call get_time (time%model_time, isec2, iday2)
daymodel = iday2 + isec2/86400.
end_of_month = set_switch (dayint, time%model_time, dt_time)

!
!       if this is the end of month, make sure that the month and
!       year pointers point to the month/year just completed, and not
!       possibly to the next month/year. This is important for indexing
!       purposes
!

if (end_of_month) then  !{

!
!       check whether we think we're in the next month and if so,
!       decrement the month and possibly year
!

  if (days .lt. 15) then  !{

!
!       if we're in the next month then we need to handle "January"
!       differently (namely go back to December of the previous year)
!

    if (months .eq. 1) then  !{

      imonth = 12
      iyear = years - 1

    else  !}{

!
!       otherwise just decrement the month
!

      imonth = months - 1
      iyear = years

    endif  !}

  else  !}{

!
!       we're think that we are at the end of the just-processed month
!       so no modifications need be done
!

    imonth = months
    iyear = years

  endif  !}

else  !}{

!
!       not end of month case, so just save the month
!

  imonth = months
  iyear = years

endif  !}

!
!       set a correct end of year indicator
!

end_of_year = end_of_month .and. imonth .eq. 12

return

contains

function set_switch (switch_interval, time_since, dt_time)  !{

!
!       Function definition
!

logical                 :: set_switch

!
!       Arguments
!

real, intent(in)                :: switch_interval      ! in units of days
type(time_type), intent(in)     :: time_since
type(time_type), intent(in)     :: dt_time

!
!       local variables
!

type(time_type)         :: interval_time
type(time_type)         :: current_time
type(time_type)         :: next_time
type(time_type)         :: half_dt_time
integer                 :: n
integer                 :: seconds
integer                 :: days

if (switch_interval < 0) then
  set_switch = .false.
  return
endif

days           = int(switch_interval)
seconds        = (switch_interval - days)*86400
interval_time  = set_time (seconds, days)

if (interval_time <= dt_time) then
  set_switch = .true.
else
  half_dt_time = dt_time / 2
  n            = (time_since + half_dt_time) / interval_time
  current_time = time_since - n * interval_time
  if (current_time  <= half_dt_time) then
    next_time = (time_since + dt_time) - n * interval_time
    if (current_time == next_time) then
      set_switch = .false.
    else
      set_switch = .true.
    endif
  else
    set_switch = .false.
  endif
endif

end function  set_switch  !}

end subroutine do_time_calc  !}
! </SUBROUTINE> NAME="do_time_calc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_bbc">
!
! <DESCRIPTION>
!       call subroutines to perform bottom boundary condition
!       calculations
! </DESCRIPTION>

subroutine ocean_tpm_bbc  !{

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_bbc'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'

!
!       set some indices and flags dependent on time
!

call do_time_calc

! JornB
if (do_ocean_tracer_fabm) then  !{
  call ocean_tracer_fabm_bbc
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_bbc
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_bbc
endif  !}

if (do_ocean_bgc_phyto) then  !{
  call ocean_bgc_phyto_bbc
endif  !}

return

end subroutine ocean_tpm_bbc  !}
! </SUBROUTINE> NAME="ocean_tpm_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_end">
!
! <DESCRIPTION>
!       Finish up calculations for the tracer packages,
!       possibly writing out non-field restart information
! </DESCRIPTION>
!

subroutine ocean_tpm_end(Thickness)  !{

type(ocean_thickness_type), intent(in) :: Thickness

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_end'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'

!
!       call subroutines to finish up the run
!

if (do_ocean_age_tracer) then  !{
  call ocean_age_tracer_end
endif  !}

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_end
endif  !}

! JornB
if (do_ocean_tracer_fabm) then  !{
  call ocean_tracer_fabm_end(Thickness)
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_end(Thickness)
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_end(Thickness)
endif  !}

if (do_ocean_bgc_phyto) then  !{
  call ocean_bgc_phyto_end(Thickness)
endif  !}

return

end subroutine ocean_tpm_end  !}
! </SUBROUTINE> NAME="ocean_tpm_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_sbc">
!
! <DESCRIPTION>
!       call subroutines to perform surface boundary condition
!       calculations
! </DESCRIPTION>
!

subroutine ocean_tpm_sbc(robert)  !{

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_sbc'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'

!
!       Arguments
!

real, intent(in)        :: robert       ! robert time-filter coefficient

!
!       set some indices and flags dependent on time
!

call do_time_calc

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_sbc
endif  !}

! JornB
if (do_ocean_tracer_fabm) then  !{
  call ocean_tracer_fabm_sbc(robert)
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_sbc(robert)
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_sbc(robert)
endif  !}

if (do_ocean_bgc_phyto) then  !{
  call ocean_bgc_phyto_sbc(robert)
endif  !}

return

end subroutine ocean_tpm_sbc  !}
! </SUBROUTINE> NAME="ocean_tpm_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>
!

subroutine ocean_tpm_init(Domain_in, Grid_in, Time_in,          &
                          dtts_in)  !{

!use ocean_tpm_pointers_mod
use ocean_tpm_util_mod, only: domain, grid, time, dtts, isc, iec, jsc, jec, nk,       &
     isd, ied, jsd, jed, indtemp, indsal
use field_manager_mod,           only: fm_get_index
use ocean_types_mod,             only: ocean_domain_type, ocean_grid_type
use ocean_types_mod,             only: ocean_time_type, ocean_prog_tracer_type
!
!-----------------------------------------------------------------------
!       arguments
!-----------------------------------------------------------------------
!

type(ocean_domain_type), intent(in), target                     :: domain_in
type(ocean_grid_type), intent(in), target                       :: grid_in
type(ocean_time_type), intent(in), target                       :: time_in
real, intent(in)                                                :: dtts_in

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_init'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!


!
!-----------------------------------------------------------------------
!       save pointers to various model types so that they may be used
!       by other tracer packages
!-----------------------------------------------------------------------
!

domain => domain_in
grid   => grid_in
time   => time_in
!t_prog => t_prog_in
!t_diag => t_diag_in
dtts   =  dtts_in

!
!       Save some indices which are used often
!

isc = domain%isc
iec = domain%iec
jsc = domain%jsc
jec = domain%jec
isd = domain%isd
ied = domain%ied
jsd = domain%jsd
jed = domain%jed
nk  = grid%nk

!
!       Determine indices for temperature and salinity
!

indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
if (indtemp .le. 0) then  !{
  call mpp_error(FATAL,trim(error_header) // ' Could not get the temperature index')
endif  !}

indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
if (indsal .le. 0) then  !{
  call mpp_error(FATAL,trim(error_header) // ' Could not get the salinity index')
endif  !}

!
!-----------------------------------------------------------------------
!       Check which tracer packages have been turned on
!-----------------------------------------------------------------------
!

!
!       Call subroutines to perform initialization operations
!

call ocean_age_tracer_init

call ocmip2_cfc_init

! JornB
call ocean_tracer_fabm_init

call ocmip2_biotic_init

call ocean_bgc_restore_init

call ocean_bgc_phyto_init

return

end subroutine ocean_tpm_init  !}
! </SUBROUTINE> NAME="ocean_tpm_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_source">
!
! <DESCRIPTION>
!       Calculate the source arrays for the tracer packages
! </DESCRIPTION>
!

subroutine ocean_tpm_source(Thickness)

 type(ocean_thickness_type), intent(in) :: Thickness
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_source'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'

!
!       set some indices and flags dependent on time
!

call do_time_calc

!
!       Call subroutines to determine the source array
!

if (do_ocean_age_tracer) then
  call ocean_age_tracer_source(Thickness)
endif

! JornB
if (do_ocean_tracer_fabm) then
  call ocean_tracer_fabm_source(Thickness)
endif

if (do_ocmip2_biotic) then
  call ocmip2_biotic_source(Thickness)
endif

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_source(Thickness)
endif  !}

if (do_ocean_bgc_phyto) then  !{
  call ocean_bgc_phyto_source(Thickness)
endif  !}

return

end subroutine ocean_tpm_source
! </SUBROUTINE> NAME="ocean_tpm_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_start">
!
! <DESCRIPTION>
!       Start the tracer packages.
!       This could include reading in extra restart information,
!       processing namelists or doing initial calculations
! </DESCRIPTION>
!

subroutine ocean_tpm_start(T_prog_in, T_diag_in,Dens)  !{

use ocean_tpm_util_mod, only: t_prog, t_diag
use ocean_types_mod, only: ocean_prog_tracer_type, ocean_diag_tracer_type

type(ocean_density_type), target, intent(in) :: Dens

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_start'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       arguments
!-----------------------------------------------------------------------
!

type(ocean_prog_tracer_type), intent(in), target, dimension(:)  :: t_prog_in
type(ocean_diag_tracer_type), intent(in), target, dimension(:)  :: t_diag_in

!
!-----------------------------------------------------------------------
!       save pointers to various model types so that they may be used
!       by other tracer packages
!-----------------------------------------------------------------------
!

t_prog => t_prog_in
t_diag => t_diag_in

!
!       call subroutines to start the tracer packages
!

if (do_ocean_age_tracer) then  !{
  call ocean_age_tracer_start
endif  !}

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_start
endif  !}

! JornB
if (do_ocean_tracer_fabm) then  !{
  call ocean_tracer_fabm_start(Dens)
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_start
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_start
endif  !}

if (do_ocean_bgc_phyto) then  !{
  call ocean_bgc_phyto_start
endif  !}

return

end subroutine ocean_tpm_start  !}
! </SUBROUTINE> NAME="ocean_tpm_start"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_tracer">
!
! <DESCRIPTION>
!       Subroutine to do calculations needed every time-step after
!       the continuity equation has been integrated
! </DESCRIPTION>
!

subroutine ocean_tpm_tracer  !{

use ocean_tpm_util_mod, only: dtts

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_tracer'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'

!
!       call subroutines to perform functions required each time-step
!       after the continuity equation has been integrated
!

if (do_ocean_age_tracer) then  !{
  call ocean_age_tracer_tracer
endif  !}

! JornB
if (do_ocean_tracer_fabm) then  !{
  call ocean_tracer_fabm_tracer
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_tracer
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_tracer
endif  !}

if (do_ocean_bgc_phyto) then  !{
  call ocean_bgc_phyto_tracer
endif  !}

return

end subroutine ocean_tpm_tracer  !}
! </SUBROUTINE> NAME="ocean_tpm_tracer"

end module ocean_tpm_mod  !}

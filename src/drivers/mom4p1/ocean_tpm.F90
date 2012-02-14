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
!       ocean_tpm_flux_init: this routine initalizes field elements
!               relating to the ocean-atmosphere gas fluxes
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
!
!       The following routines are called in relation to tying in to
!               the FMS coupler to calculate fluxes for the additional
!               tracers:
!
!       ocean_tpm_init_sfc: Allocate arrays for the accumulation of
!               data to be used by the coupler
!
!       ocean_tpm_sum_sfc: Accumulate data for the coupler
!
!       ocean_tpm_avg_sfc: Take the time-mean of the fields for the coupler
!
!       ocean_tpm_zero_sfc: Zero out the fields for the coupler to allow
!               for accumulation for the next time period
!
!       ocean_tpm_sfc_end: Save out fields for the restart.
!
!</DESCRIPTION>
!
! <INFO>
! </INFO>
!

use mpp_mod,            only: stdout, mpp_error, FATAL
use mpp_domains_mod,    only: mpp_get_compute_domain
use ocean_types_mod,    only: ocean_thickness_type, ocean_public_type
use ocean_types_mod,    only: ocean_options_type
use ocean_types_mod,    only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,    only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,    only: ocean_density_type

!
!       Place tracer modules here
!

use ocean_age_tracer_mod, only: do_ocean_age_tracer
use ocean_age_tracer_mod, only: ocean_age_tracer_init
use ocean_age_tracer_mod, only: ocean_age_tracer_source
use ocean_age_tracer_mod, only: ocean_age_tracer_start
use ocean_age_tracer_mod, only: ocean_age_tracer_tracer

use ocean_residency_mod, only: do_ocean_residency
use ocean_residency_mod, only: ocean_residency_init
use ocean_residency_mod, only: ocean_residency_source
use ocean_residency_mod, only: ocean_residency_start
use ocean_residency_mod, only: ocean_residency_tracer

#ifdef USE_OCEAN_BGC
use ocean_pert_co2_mod, only: do_ocean_pert_co2
use ocean_pert_co2_mod, only: ocean_pert_co2_avg_sfc
use ocean_pert_co2_mod, only: ocean_pert_co2_end
use ocean_pert_co2_mod, only: ocean_pert_co2_flux_init
use ocean_pert_co2_mod, only: ocean_pert_co2_init
use ocean_pert_co2_mod, only: ocean_pert_co2_init_sfc
use ocean_pert_co2_mod, only: ocean_pert_co2_sbc
use ocean_pert_co2_mod, only: ocean_pert_co2_source
use ocean_pert_co2_mod, only: ocean_pert_co2_sum_sfc
use ocean_pert_co2_mod, only: ocean_pert_co2_start
use ocean_pert_co2_mod, only: ocean_pert_co2_zero_sfc

use ocmip2_abiotic_mod, only: do_ocmip2_abiotic
use ocmip2_abiotic_mod, only: ocmip2_abiotic_avg_sfc
use ocmip2_abiotic_mod, only: ocmip2_abiotic_end
use ocmip2_abiotic_mod, only: ocmip2_abiotic_flux_init
use ocmip2_abiotic_mod, only: ocmip2_abiotic_init
use ocmip2_abiotic_mod, only: ocmip2_abiotic_init_sfc
use ocmip2_abiotic_mod, only: ocmip2_abiotic_sbc
use ocmip2_abiotic_mod, only: ocmip2_abiotic_sfc_end
use ocmip2_abiotic_mod, only: ocmip2_abiotic_source
use ocmip2_abiotic_mod, only: ocmip2_abiotic_sum_sfc
use ocmip2_abiotic_mod, only: ocmip2_abiotic_start
use ocmip2_abiotic_mod, only: ocmip2_abiotic_zero_sfc
use ocmip2_abiotic_mod, only: ocmip2_abiotic_restart
use ocmip2_abiotic_mod, only: ocmip2_abiotic_tracer

use ocmip2_cfc_mod, only: do_ocmip2_cfc
use ocmip2_cfc_mod, only: ocmip2_cfc_avg_sfc
use ocmip2_cfc_mod, only: ocmip2_cfc_end
use ocmip2_cfc_mod, only: ocmip2_cfc_flux_init
use ocmip2_cfc_mod, only: ocmip2_cfc_init
use ocmip2_cfc_mod, only: ocmip2_cfc_init_sfc
use ocmip2_cfc_mod, only: ocmip2_cfc_sbc
use ocmip2_cfc_mod, only: ocmip2_cfc_sfc_end
use ocmip2_cfc_mod, only: ocmip2_cfc_start
use ocmip2_cfc_mod, only: ocmip2_cfc_sum_sfc
use ocmip2_cfc_mod, only: ocmip2_cfc_zero_sfc

use ocmip2_biotic_mod, only: do_ocmip2_biotic
use ocmip2_biotic_mod, only: ocmip2_biotic_avg_sfc
use ocmip2_biotic_mod, only: ocmip2_biotic_bbc
use ocmip2_biotic_mod, only: ocmip2_biotic_end
use ocmip2_biotic_mod, only: ocmip2_biotic_flux_init
use ocmip2_biotic_mod, only: ocmip2_biotic_init
use ocmip2_biotic_mod, only: ocmip2_biotic_init_sfc
use ocmip2_biotic_mod, only: ocmip2_biotic_sbc
use ocmip2_biotic_mod, only: ocmip2_biotic_sfc_end
use ocmip2_biotic_mod, only: ocmip2_biotic_source
use ocmip2_biotic_mod, only: ocmip2_biotic_sum_sfc
use ocmip2_biotic_mod, only: ocmip2_biotic_start
use ocmip2_biotic_mod, only: ocmip2_biotic_zero_sfc
use ocmip2_biotic_mod, only: ocmip2_biotic_restart

use ocean_bgc_restore_mod, only: do_ocean_bgc_restore
use ocean_bgc_restore_mod, only: ocean_bgc_restore_avg_sfc
use ocean_bgc_restore_mod, only: ocean_bgc_restore_bbc
use ocean_bgc_restore_mod, only: ocean_bgc_restore_end
use ocean_bgc_restore_mod, only: ocean_bgc_restore_flux_init
use ocean_bgc_restore_mod, only: ocean_bgc_restore_init
use ocean_bgc_restore_mod, only: ocean_bgc_restore_init_sfc
use ocean_bgc_restore_mod, only: ocean_bgc_restore_sbc
use ocean_bgc_restore_mod, only: ocean_bgc_restore_sfc_end
use ocean_bgc_restore_mod, only: ocean_bgc_restore_source
use ocean_bgc_restore_mod, only: ocean_bgc_restore_sum_sfc
use ocean_bgc_restore_mod, only: ocean_bgc_restore_start
use ocean_bgc_restore_mod, only: ocean_bgc_restore_zero_sfc
use ocean_bgc_restore_mod, only: ocean_bgc_restore_restart

use ocmip2_he_mod, only: do_ocmip2_he
use ocmip2_he_mod, only: ocmip2_he_avg_sfc
use ocmip2_he_mod, only: ocmip2_he_end
use ocmip2_he_mod, only: ocmip2_he_flux_init
use ocmip2_he_mod, only: ocmip2_he_init
use ocmip2_he_mod, only: ocmip2_he_init_sfc
use ocmip2_he_mod, only: ocmip2_he_sbc
use ocmip2_he_mod, only: ocmip2_he_source
use ocmip2_he_mod, only: ocmip2_he_start
use ocmip2_he_mod, only: ocmip2_he_sum_sfc
use ocmip2_he_mod, only: ocmip2_he_zero_sfc
use ocmip2_he_mod, only: ocmip2_he_restart

use ocean_po4_pre_mod, only: do_ocean_po4_pre
use ocean_po4_pre_mod, only: ocean_po4_pre_end
use ocean_po4_pre_mod, only: ocean_po4_pre_init
use ocean_po4_pre_mod, only: ocean_po4_pre_start
use ocean_po4_pre_mod, only: ocean_po4_pre_tracer
use ocean_po4_pre_mod, only: ocean_po4_pre_zero_sfc

use ocean_ibgc_mod, only: do_ocean_ibgc
use ocean_ibgc_mod, only: ocean_ibgc_avg_sfc
use ocean_ibgc_mod, only: ocean_ibgc_bbc
use ocean_ibgc_mod, only: ocean_ibgc_end
use ocean_ibgc_mod, only: ocean_ibgc_flux_init
use ocean_ibgc_mod, only: ocean_ibgc_init
use ocean_ibgc_mod, only: ocean_ibgc_init_sfc
use ocean_ibgc_mod, only: ocean_ibgc_sbc
use ocean_ibgc_mod, only: ocean_ibgc_sfc_end
use ocean_ibgc_mod, only: ocean_ibgc_source
use ocean_ibgc_mod, only: ocean_ibgc_sum_sfc
use ocean_ibgc_mod, only: ocean_ibgc_start
use ocean_ibgc_mod, only: ocean_ibgc_tracer
use ocean_ibgc_mod, only: ocean_ibgc_zero_sfc
use ocean_ibgc_mod, only: ocean_ibgc_restart

use ocean_generic_mod, only: do_generic_tracer
use ocean_generic_mod, only: ocean_generic_sum_sfc
use ocean_generic_mod, only: ocean_generic_zero_sfc
use ocean_generic_mod, only: ocean_generic_sbc
use ocean_generic_mod, only: ocean_generic_init
use ocean_generic_mod, only: ocean_generic_column_physics
use ocean_generic_mod, only: ocean_generic_end
use ocean_generic_mod, only: ocean_generic_flux_init

use ocean_fabm_mod, only: do_ocean_fabm
use ocean_fabm_mod, only: ocean_fabm_end
use ocean_fabm_mod, only: ocean_fabm_init
use ocean_fabm_mod, only: ocean_fabm_source
use ocean_fabm_mod, only: ocean_fabm_start
use ocean_fabm_mod, only: ocean_fabm_bbc
use ocean_fabm_mod, only: ocean_fabm_sbc
use ocean_fabm_mod, only: ocean_fabm_tracer
use ocean_fabm_mod, only: ocean_fabm_avg_sfc
use ocean_fabm_mod, only: ocean_fabm_sum_sfc
use ocean_fabm_mod, only: ocean_fabm_zero_sfc
use ocean_fabm_mod, only: ocean_fabm_flux_init
use ocean_fabm_mod, only: ocean_fabm_sfc_end
use ocean_fabm_mod, only: ocean_fabm_init_sfc

#endif

use ocean_frazil_mod,     only: ocean_frazil_init
use ocean_tempsalt_mod,   only: ocean_tempsalt_init
use ocean_passive_mod,    only: ocean_passive_init

use transport_matrix_mod, only: do_transport_matrix
use transport_matrix_mod, only: transport_matrix_init
use transport_matrix_mod, only: transport_matrix_start
use transport_matrix_mod, only: transport_matrix_store_implicit

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
public ocean_tpm_flux_init
public ocean_tpm_sbc
public ocean_tpm_source
public ocean_tpm_start
public ocean_tpm_tracer
public ocean_tpm_init_sfc
public ocean_tpm_sum_sfc
public ocean_tpm_avg_sfc
public ocean_tpm_zero_sfc
public ocean_tpm_sfc_end

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

character(len=128) :: version = '$Id: ocean_tpm.F90,v 16.0.2.2.32.2.12.7 2009/10/14 14:37:24 smg Exp $'
character(len=128) :: tagname = '$Name: mom4p1_pubrel_dec2009_nnz $'

integer :: imonth
integer :: iyear
logical :: end_of_day
logical :: end_of_month
logical :: end_of_year
logical :: mid_month

contains

!#######################################################################
! <SUBROUTINE NAME="do_time_calc">
!
! <DESCRIPTION>
!       call subroutines to perform time calculations
! </DESCRIPTION>
!

subroutine do_time_calc(time, dtts)  !{

!
!-----------------------------------------------------------------------
!
!       Modules
!
!-----------------------------------------------------------------------
!

use time_manager_mod, only: time_type, set_time, get_time
use time_manager_mod, only: set_date, get_date, days_in_month, operator(+)
use time_manager_mod, only: operator(<=), operator(==)
use time_manager_mod, only: operator(*), operator(/), operator(-)

implicit none

!
!-----------------------------------------------------------------------
!     Arguments
!-----------------------------------------------------------------------
!

type(ocean_time_type), intent(in)       :: time
real, intent(in)                        :: dtts

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
! Check that old ifdef is not accidently defined
#ifdef USE_OCEAN_OCMIP2
  call mpp_error(FATAL, &
  '==>Error in ocean_tmp_mod: cpp option USE_OCEAN_OCMIP2 is now called USE_OCEAN_BGC. Please recompile...sorry.')
#endif


time_tau = time%tau

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

implicit none

!
!       Function definition
!

logical                 :: set_switch

!
!----------------------------------------------------------------------
!       Arguments
!----------------------------------------------------------------------
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

subroutine ocean_tpm_bbc(Domain, Grid, T_prog)  !{

!
!-----------------------------------------------------------------------
!     Arguments
!-----------------------------------------------------------------------
!

type(ocean_domain_type), intent(in)                             :: Domain
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog

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


#ifdef USE_OCEAN_BGC

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_bbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec,        &
                         Domain%isd, Domain%ied, Domain%jsd, Domain%jed, T_prog, Grid%kmt)
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_bbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec,    &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed, T_prog, Grid%kmt)
endif  !}

if (do_ocean_fabm) then  !{
  call ocean_fabm_bbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec,        &
                         Domain%isd, Domain%ied, Domain%jsd, Domain%jed, T_prog, Grid%kmt)
endif  !}

#endif

return

end subroutine ocean_tpm_bbc  !}
! </SUBROUTINE> NAME="ocean_tpm_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ocean_tpm_restart(time_stamp)
  character(len=*),           intent(in), optional :: time_stamp

#ifdef USE_OCEAN_BGC

if (do_ocmip2_abiotic) call ocmip2_abiotic_restart(time_stamp)
if (do_ocmip2_biotic) call ocmip2_biotic_restart(time_stamp)
if (do_ocmip2_he) call ocmip2_he_restart(time_stamp)
if (do_ocean_bgc_restore) call ocean_bgc_restore_restart(time_stamp)
if (do_ocean_ibgc) call ocean_ibgc_restart(time_stamp)
!if (do_ocean_fabm) call ocean_fabm_restart(time_stamp)

#endif

end subroutine ocean_tpm_restart
! </SUBROUTINE> NAME="ocean_tpm_restart"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_end">
!
! <DESCRIPTION>
!       Finish up calculations for the tracer packages,
!       possibly writing out non-field restart information
! </DESCRIPTION>
!

subroutine ocean_tpm_end(Domain, Grid, T_prog, T_diag, Time, Thickness)  !{

implicit none

!
!-----------------------------------------------------------------------
!     Arguments
!-----------------------------------------------------------------------
!

type(ocean_domain_type), intent(in)                     :: Domain
type(ocean_grid_type), intent(in)                       :: Grid
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
type(ocean_diag_tracer_type), dimension(:), intent(in)  :: T_diag
type(ocean_time_type), intent(in)                       :: Time
type(ocean_thickness_type), intent(in)                  :: Thickness

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

#ifdef USE_OCEAN_BGC

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_end(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,  &
                      Domain%isd, Domain%ied, Domain%jsd, Domain%jed,           &
                      T_prog, Grid%dat, Grid%tmask,                             &
                      Thickness%rho_dzt, Time%taup1)
endif  !}

if (do_ocean_pert_co2) then  !{
  call ocean_pert_co2_end(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,      &
                          Domain%isd, Domain%ied, Domain%jsd, Domain%jed,               &
                          T_prog, grid%dat, grid%tmask, Domain%domain2d,                &
                          Thickness%rho_dzt, Time%taup1)
endif  !}

if (do_ocmip2_abiotic) then  !{
  call ocmip2_abiotic_end(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,      &
                          Domain%isd, Domain%ied, Domain%jsd, Domain%jed,               &
                          T_prog, Grid%dat, Grid%tmask, Domain%domain2d,                &
                          Thickness%rho_dzt, Time%taup1)
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_end(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,       &
                         Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                &
                         T_prog, Grid%dat, Grid%tmask, Domain%domain2d,                 &
                         Thickness%rho_dzt, Time%taup1)
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_end(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,   &
                           Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                           T_prog, T_diag, Grid%dat, Grid%tmask, Domain%domain2d,       &
                           Thickness%rho_dzt, Time%taup1)
endif  !}

if (do_ocmip2_he) then  !{
  call ocmip2_he_end(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,  &
                     Domain%isd, Domain%ied, Domain%jsd, Domain%jed,           &
                     T_prog, grid%dat, grid%tmask,                             &
                     Thickness%rho_dzt, Time%taup1)
endif  !}

if (do_ocean_po4_pre) then  !{
  call ocean_po4_pre_end(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk, &
                     Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                     T_prog, Grid%dat, Grid%tmask, Thickness%rho_dzt, Time%taup1)
endif  !}

if (do_ocean_ibgc) then  !{
  call ocean_ibgc_end(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk, &
                     Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                     T_prog, T_diag, Grid%dat, Grid%tmask, Domain%domain2d,        &
                     Thickness%rho_dzt, Time%taup1)
endif  !}

if (do_generic_tracer) call ocean_generic_end

if (do_ocean_fabm) then  !{
  call ocean_fabm_end(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,       &
                         Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                &
                         T_prog, Grid%dat, Grid%tmask, Domain%domain2d,                 &
                         Thickness%rho_dzt, Time%taup1)
endif  !}

#endif

return

end subroutine ocean_tpm_end  !}
! </SUBROUTINE> NAME="ocean_tpm_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_init_sfc">
!
! <DESCRIPTION>
!       call subroutines to perform surface coupler initializations
!
!       Note: this subroutine should be merged into ocean_tpm_start
! </DESCRIPTION>
!

subroutine ocean_tpm_init_sfc(Domain, T_prog, Dens, Ocean, Time, Grid)  !{

implicit none

!
!-----------------------------------------------------------------------
!     Arguments
!-----------------------------------------------------------------------
!

type(ocean_domain_type), intent(in)                     :: Domain
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
type(ocean_density_type), intent(in)                    :: Dens
type(ocean_public_type), intent(inout)                  :: Ocean
type(ocean_time_type), intent(in)                       :: Time
type(ocean_grid_type), intent(in)                       :: Grid

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_init_sfc'
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

logical, save   :: initialized = .false.
integer         :: isc_bnd
integer         :: iec_bnd
integer         :: jsc_bnd
integer         :: jec_bnd

if (.not. initialized) then  !{

  call mpp_get_compute_domain(Ocean%Domain, isc_bnd, iec_bnd, jsc_bnd, jec_bnd)

#ifdef USE_OCEAN_BGC

  if (do_ocmip2_cfc) then  !{
    call ocmip2_cfc_init_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,           &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                    &
                             isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                &
                             Ocean%fields, T_prog, Dens%rho, Time%taum1,                        &
                             Grid%tmask)
  endif  !}

  if (do_ocean_pert_co2) then  !{
    call ocean_pert_co2_init_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,       &
                                 Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                &
                                 isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                            &
                                 Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,   &
                                 Grid%tmask)
  endif  !}

  if (do_ocmip2_abiotic) then  !{
    call ocmip2_abiotic_init_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,       &
                                 Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                &
                                 isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                            &
                                 Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,   &
                                 Grid%tmask)
  endif  !}

  if (do_ocmip2_biotic) then  !{
    call ocmip2_biotic_init_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,        &
                                Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                 &
                                isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                             &
                                Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,    &
                                Grid%tmask)
  endif  !}

  if (do_ocean_bgc_restore) then  !{
    call ocean_bgc_restore_init_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,    &
                                    Domain%isd, Domain%ied, Domain%jsd, Domain%jed,             &
                                    isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                         &
                                    Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,&
                                    Grid%tmask)
  endif  !}

  if (do_generic_tracer) call ocean_generic_sum_sfc(Domain%isd,Domain%jsd, Ocean, T_prog, Dens, Time )

  if (do_ocmip2_he) then  !{
    call ocmip2_he_init_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,   &
                            Domain%isd, Domain%ied, Domain%jsd, Domain%jed,            &
                            isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                        &
                            Ocean%fields, T_prog, Dens%rho, Time%taum1, Grid%tmask)
  endif  !}

  if (do_ocean_ibgc) then  !{
    call ocean_ibgc_init_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,      &
                              Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                   &
                              isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                               &
                              Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,      &
                              Grid%tmask)
  endif  !}

  if (do_ocean_fabm) then  !{
    call ocean_fabm_init_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,      &
                              Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                   &
                              isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                               &
                              Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,      &
                              Grid%tmask)
  endif  !}

#endif


  initialized = .true.

endif  !}

return

end subroutine ocean_tpm_init_sfc  !}
! </SUBROUTINE> NAME="ocean_tpm_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_sum_sfc">
!
! <DESCRIPTION>
!       call subroutines to perform surface coupler initializations
! </DESCRIPTION>
!

subroutine ocean_tpm_sum_sfc(Domain, T_prog, Dens, Ocean, Time, Grid, isc_bnd, iec_bnd, jsc_bnd, jec_bnd)  !{

implicit none

!
!-----------------------------------------------------------------------
!     Arguments
!-----------------------------------------------------------------------
!

type(ocean_domain_type), intent(in)                     :: Domain
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
type(ocean_density_type), intent(in)                    :: Dens
type(ocean_public_type), intent(inout)                  :: Ocean
type(ocean_time_type), intent(in)                       :: Time
type(ocean_grid_type), intent(in)                       :: Grid
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: iec_bnd
integer, intent(in)                                     :: jsc_bnd
integer, intent(in)                                     :: jec_bnd

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_sum_sfc'
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

#ifdef USE_OCEAN_BGC

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_sum_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,            &
                          Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                     &
                          isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                 &
                          Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,        &
                          Grid%tmask)
endif  !}

if (do_ocean_pert_co2) then  !{
  call ocean_pert_co2_sum_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,          &
                              Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                   &
                              isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                               &
                              Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,      &
                              Grid%tmask)
endif  !}

if (do_ocmip2_abiotic) then  !{
  call ocmip2_abiotic_sum_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,        &
                              Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                 &
                              isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                             &
                              Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,    &
                              Grid%tmask)
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_sum_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,         &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                  &
                             isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                              &
                             Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,     &
                             Grid%tmask)
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_sum_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,     &
                                 Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                                 isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                          &
                                 Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time, &
                                 Grid%tmask)
endif  !}

if (do_generic_tracer) then
   call ocean_generic_sum_sfc(Domain%isd,Domain%jsd, Ocean, T_prog, Dens, Time )
endif

if (do_ocmip2_he) then  !{
  call ocmip2_he_sum_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,             &
                          Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                     &
                          isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                 &
                          Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,        &
                          Grid%tmask)
endif  !}

if (do_ocean_ibgc) then  !{
  call ocean_ibgc_sum_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,       &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                  &
                             isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                              &
                             Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,     &
                             Grid%tmask)
endif  !}

if (do_ocean_fabm) then  !{
  call ocean_fabm_sum_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,       &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                  &
                             isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                              &
                             Ocean%fields, T_prog, Dens%rho, Time%taum1, Time%model_time,     &
                             Grid%tmask)
endif  !}

#endif

return

end subroutine ocean_tpm_sum_sfc  !}
! </SUBROUTINE> NAME="ocean_tpm_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_avg_sfc">
!
! <DESCRIPTION>
!       call subroutines to perform surface coupler initializations
! </DESCRIPTION>
!

subroutine ocean_tpm_avg_sfc(Domain, Ocean, Grid, isc_bnd, iec_bnd, jsc_bnd, jec_bnd)  !{

implicit none

!
!-----------------------------------------------------------------------
!     Arguments
!-----------------------------------------------------------------------
!

type(ocean_domain_type), intent(in)     :: Domain
type(ocean_public_type), intent(inout)  :: Ocean
type(ocean_grid_type), intent(in)       :: Grid
integer, intent(in)                     :: isc_bnd
integer, intent(in)                     :: iec_bnd
integer, intent(in)                     :: jsc_bnd
integer, intent(in)                     :: jec_bnd

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_avg_sfc'
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

#ifdef USE_OCEAN_BGC

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_avg_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,        &
                          Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                 &
                          isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                             &
                          Ocean%fields, Ocean%avg_kount, Grid%tmask)
endif  !}

if (do_ocean_pert_co2) then  !{
  call ocean_pert_co2_avg_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,  &
                              Domain%isd, Domain%ied, Domain%jsd, Domain%jed,           &
                              isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                       &
                              Ocean%fields, Ocean%avg_kount, Grid%tmask)
endif  !}

if (do_ocmip2_abiotic) then  !{
  call ocmip2_abiotic_avg_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,    &
                              Domain%isd, Domain%ied, Domain%jsd, Domain%jed,             &
                              isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                         &
                              Ocean%fields, Ocean%avg_kount, Grid%tmask)
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_avg_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,     &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                             isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                          &
                             Ocean%fields, Ocean%avg_kount, Grid%tmask)
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_avg_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk, &
                                 Domain%isd, Domain%ied, Domain%jsd, Domain%jed,          &
                                 isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                      &
                                 Ocean%fields, Ocean%avg_kount, Grid%tmask)
endif  !}

if (do_ocmip2_he) then  !{
  call ocmip2_he_avg_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,         &
                          Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                 &
                          isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                             &
                          Ocean%fields, Ocean%avg_kount, Grid%tmask)
endif  !}

if (do_ocean_ibgc) then  !{
  call ocean_ibgc_avg_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,   &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                             isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                          &
                             Ocean%fields, Ocean%avg_kount, Grid%tmask)
endif  !}

if (do_ocean_fabm) then  !{
  call ocean_fabm_avg_sfc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,   &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                             isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                          &
                             Ocean%fields, Ocean%avg_kount, Grid%tmask)
endif  !}


#endif


return

end subroutine ocean_tpm_avg_sfc  !}
! </SUBROUTINE> NAME="ocean_tpm_avg_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_zero_sfc">
!
! <DESCRIPTION>
!       call subroutines to perform surface coupler initializations
! </DESCRIPTION>
!

subroutine ocean_tpm_zero_sfc(Ocean)  !{

implicit none

!
!-----------------------------------------------------------------------
!     Arguments
!-----------------------------------------------------------------------
!

type(ocean_public_type), intent(inout)  :: Ocean

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_zero_sfc'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'


#ifdef USE_OCEAN_BGC

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_zero_sfc(Ocean%fields)
endif  !}

if (do_ocmip2_he) then  !{
  call ocmip2_he_zero_sfc(Ocean%fields)
endif  !}

if (do_ocean_pert_co2) then  !{
  call ocean_pert_co2_zero_sfc(Ocean%fields)
endif  !}

if (do_ocmip2_abiotic) then  !{
  call ocmip2_abiotic_zero_sfc(Ocean%fields)
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_zero_sfc(Ocean%fields)
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_zero_sfc(Ocean%fields)
endif  !}

if (do_ocean_ibgc) then  !{
  call ocean_ibgc_zero_sfc(Ocean%fields)
endif  !}

if (do_generic_tracer) call ocean_generic_zero_sfc(Ocean%fields)

if (do_ocean_fabm) then  !{
  call ocean_fabm_zero_sfc(Ocean%fields)
endif  !}

#endif

return

end subroutine ocean_tpm_zero_sfc  !}
! </SUBROUTINE> NAME="ocean_tpm_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_sfc_end">
!
! <DESCRIPTION>
!       call subroutines to perform surface coupler initializations
! </DESCRIPTION>
!

subroutine ocean_tpm_sfc_end !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_sfc_end'
character(len=256), parameter   :: error_header = '==>Error from ' // trim(mod_name) //   &
                                                  '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(mod_name) //  &
                                                 '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header = '==>Note from ' // trim(mod_name) //     &
                                                 '(' // trim(sub_name) // '):'


#ifdef USE_OCEAN_BGC

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_sfc_end
endif  !}

if (do_ocmip2_abiotic) then  !{
  call ocmip2_abiotic_sfc_end
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_sfc_end
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_sfc_end
endif  !}

if (do_ocean_ibgc) then  !{
  call ocean_ibgc_sfc_end
endif  !}

if (do_ocean_fabm) then  !{
  call ocean_fabm_sfc_end
endif  !}

#endif


return

end subroutine ocean_tpm_sfc_end  !}
! </SUBROUTINE> NAME="ocean_tpm_sfc_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_sbc">
!
! <DESCRIPTION>
!       call subroutines to perform surface boundary condition
!       calculations
! </DESCRIPTION>
!

subroutine ocean_tpm_sbc(Domain, Grid, T_prog, Time, Ice_ocean_boundary_fluxes, &
     runoff, isc_bnd, iec_bnd, jsc_bnd, jec_bnd)  !{

use coupler_types_mod, only: coupler_2d_bc_type

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

type(ocean_domain_type), intent(in)                             :: Domain
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_time_type), intent(in)                               :: Time
type(coupler_2d_bc_type), intent(in)                            :: Ice_ocean_boundary_fluxes
integer, intent(in)                                             :: isc_bnd
integer, intent(in)                                             :: iec_bnd
integer, intent(in)                                             :: jsc_bnd
integer, intent(in)                                             :: jec_bnd
real, dimension(Domain%isd:,Domain%jsd:), intent(in)            :: runoff

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
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!


#ifdef USE_OCEAN_BGC

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_sbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,          &
                      Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                   &
                      isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                               &
                      T_prog, Time%model_time,                                          &
                      Grid%tmask, Ice_ocean_boundary_fluxes)
endif  !}

if (do_ocean_pert_co2) then  !{
  call ocean_pert_co2_sbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,      &
                          Domain%isd, Domain%ied, Domain%jsd, Domain%jed,               &
                          isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                           &
                          T_prog, Time%taum1, Time%model_time,                          &
                          grid%tmask, Ice_ocean_boundary_fluxes)
endif  !}

if (do_ocmip2_abiotic) then  !{
  call ocmip2_abiotic_sbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,      &
                          Domain%isd, Domain%ied, Domain%jsd, Domain%jed,               &
                          isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                           &
                          T_prog, Time%taum1, Time%model_time,                          &
                          Grid%tmask, Ice_ocean_boundary_fluxes)
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_sbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,       &
                         Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                &
                         isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                            &
                         T_prog, Time%taum1, Time%model_time,                           &
                         Grid%tmask, Ice_ocean_boundary_fluxes)
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_sbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec,            &
                             Grid%nk, Domain%isd, Domain%ied, Domain%jsd,               &
                             Domain%jed, isc_bnd, iec_bnd, jsc_bnd, jec_bnd,            &
                             T_prog, Time%tau, Time%model_time,                         &
                             Grid%tmask, Ice_ocean_boundary_fluxes)
endif  !}

if (do_generic_tracer) call ocean_generic_sbc(Ice_ocean_boundary_fluxes,Domain%isd,Domain%jsd, T_prog, runoff)
if (do_ocmip2_he) then  !{
  call ocmip2_he_sbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,           &
                      Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                   &
                      isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                               &
                      T_prog, Time%model_time,                                          &
                      Grid%tmask, Ice_ocean_boundary_fluxes)
endif  !}

if (do_ocean_ibgc) then  !{
  call ocean_ibgc_sbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,     &
                           Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                           isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                          &
                           T_prog, Time%tau, Time%model_time,                           &
                           Grid%tmask, Ice_ocean_boundary_fluxes)
endif  !}

if (do_ocean_fabm) then  !{
  call ocean_fabm_sbc(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,     &
                           Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                           isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                          &
                           T_prog, Time%tau, Time%model_time,                           &
                           Grid%tmask, Ice_ocean_boundary_fluxes)
endif  !}


#endif

return

end subroutine ocean_tpm_sbc  !}
! </SUBROUTINE> NAME="ocean_tpm_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Ocean, Grid and Domains.
! </DESCRIPTION>
!

subroutine ocean_tpm_init(Domain, Grid, Time, Time_steps, &
                          Ocean_options, debug)

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

  type(ocean_domain_type),     intent(in)     :: Domain
  type(ocean_grid_type),       intent(in)     :: Grid
  type(ocean_time_type),       intent(in)     :: Time
  type(ocean_time_steps_type), intent(in)     :: Time_steps
  type(ocean_options_type),    intent(inout)  :: Ocean_options
  logical, intent(in), optional               :: debug

integer :: index_temp=-1
integer :: index_salt=-1

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
!       Check which tracer packages have been turned on
!-----------------------------------------------------------------------
!

!
!       Call subroutines to perform initialization operations
!

call ocean_tempsalt_init (Domain, Grid, Time, Ocean_options, index_temp, index_salt, debug)

call ocean_frazil_init (Domain, Grid, Time, Time_steps, Ocean_options, &
                        index_temp, index_salt, debug)

call ocean_passive_init (Domain, Grid, Ocean_options, debug)

call ocean_residency_init       ! must come first

call ocean_age_tracer_init

#ifdef USE_OCEAN_BGC

call ocmip2_cfc_init

call ocmip2_he_init

call ocean_pert_co2_init

call ocmip2_abiotic_init

call ocmip2_biotic_init

call ocean_bgc_restore_init

call ocean_po4_pre_init

call ocean_ibgc_init

call ocean_generic_init(Domain,Grid,Time)

call ocean_fabm_init

#endif

call transport_matrix_init


return

end subroutine ocean_tpm_init  !}
! </SUBROUTINE> NAME="ocean_tpm_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_tpm_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>
!

subroutine ocean_tpm_flux_init  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tpm_flux_init'
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
!       Initialize fields for the ocean-atmosphere gas fluxes
!-----------------------------------------------------------------------
!

!
!       These routines must always be called as this routine
!       is called from both atmosphere and oceanic processors to set
!       up fluxes. The variable being "fluxed" may only exist in the
!       ocean model (eg., O2), hence the need to have the calls not
!       have tests for whether the specific package is turned on--those
!       tests will always be false on atmospheric processors, even if
!       we're using that package on the oceanic processors.
!

#ifdef USE_OCEAN_BGC

call ocmip2_cfc_flux_init

call ocmip2_he_flux_init

call ocean_pert_co2_flux_init

call ocmip2_abiotic_flux_init

call ocmip2_biotic_flux_init

call ocean_bgc_restore_flux_init

call ocean_ibgc_flux_init

call ocean_generic_flux_init

call ocean_fabm_flux_init

#endif

return

end subroutine ocean_tpm_flux_init  !}
! </SUBROUTINE> NAME="ocean_tpm_flux_init"


!#######################################################################

! <SUBROUTINE NAME="ocean_tpm_source">
!
! <DESCRIPTION>
!       Calculate the source arrays for the tracer packages
! </DESCRIPTION>
!

subroutine ocean_tpm_source(isd, ied, jsd, jed, Domain, Grid, T_prog, T_diag,   &
     Time, Thickness, Dens, opacity, hblt_depth, dtts)

implicit none

!
!-----------------------------------------------------------------------
!     Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
type(ocean_domain_type), intent(in)                             :: Domain
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_diag_tracer_type), dimension(:), intent(inout)       :: T_diag
type(ocean_time_type), intent(in)                               :: Time
type(ocean_thickness_type), intent(in)                          :: Thickness
type(ocean_density_type), intent(in)                            :: Dens
real, intent(in), dimension(isd:,jsd:)                          :: hblt_depth
real, intent(in), dimension(isd:,jsd:,:)                        :: opacity
real, intent(in)                                                :: dtts
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
!       Call subroutines to determine the source array
!

if (do_ocean_age_tracer) then
  call ocean_age_tracer_source(isd, ied, jsd, jed, Grid%nk,     &
                               Time%model_time, Grid%tmask, T_prog)
endif

#ifdef USE_OCEAN_BGC
if (do_ocean_pert_co2) then
  call ocean_pert_co2_source(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,   &
                             isd, ied, jsd, jed, T_prog,                                &
                             Time%model_time, grid%tmask )
endif

if (do_ocmip2_abiotic) then
  call ocmip2_abiotic_source(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,   &
                             isd, ied, jsd, jed, T_prog,                                &
                             Time%taum1, Time%model_time, Grid%tmask,                   &
                             Thickness%rho_dzt)
endif

if (do_ocmip2_biotic) then
  call ocmip2_biotic_source(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,    &
                            isd, ied, jsd, jed, T_prog,                                 &
                            Time%taum1, Time%model_time, Grid%zw, Grid%ht, Grid%tmask,  &
                            Thickness%rho_dzt)
endif

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_source(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,&
                                isd, ied, jsd, jed, T_prog, T_diag,                     &
                                Time%taum1, Time%model_time, Grid%tmask, Grid%kmt,      &
                                Thickness%rho_dzt, dtts)
endif  !}

if (do_ocean_ibgc) then  !{
  call ocean_ibgc_source(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,  &
                              isd, ied, jsd, jed, T_prog, T_diag,                       &
                              Time%taum1, Time%model_time, Grid%dat, Grid%tmask,        &
                              Grid%kmt, Thickness%depth_zt,Dens%rho, Thickness%rho_dzt, &
                              Thickness%dzt, hblt_depth, dtts)
endif  !}

if (do_ocmip2_he) then
  call ocmip2_he_source(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,               &
                        isd, ied, jsd, jed, T_prog, Thickness%depth_zt, Thickness%dzt,         &
                        Time%taum1, Time%model_time, Grid%tmask, Grid%kmt,  Thickness%rho_dzt)
endif

if (do_ocean_fabm) then
  call ocean_fabm_source(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,    &
                            isd, ied, jsd, jed, T_prog, T_diag,                           &
                            Time%taum1, Time%model_time, Grid%dat, Grid%zw, Grid%ht, Grid%tmask,  &
                            Grid%kmt, Thickness%rho_dzt, Thickness%dzt, Dens, dtts, Domain%domain2d)
endif

#endif

if (do_ocean_residency) then  !{        ! must come last
  call ocean_residency_source(Domain%isc, Domain%iec, Domain%jsc, Domain%jec,           &
       isd, ied, jsd, jed, Grid%nk, T_prog, T_diag, Time, Thickness, Dens,              &
       grid%xt, grid%yt, grid%zw, grid%tmask, grid%kmt, hblt_depth)
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

subroutine ocean_tpm_start(Domain, Grid, T_prog, T_diag, Time, Thickness)  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

type(ocean_domain_type), intent(in)                     :: Domain
type(ocean_grid_type), intent(in)                       :: Grid
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
type(ocean_diag_tracer_type), dimension(:), intent(in)  :: T_diag
type(ocean_time_type), intent(in)                       :: Time
type(ocean_thickness_type), intent(in)                  :: Thickness

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
!       call subroutines to start the tracer packages
!

if (do_ocean_residency) then  !{        ! must come first
  call ocean_residency_start(Domain%isd, Domain%ied, Domain%jsd, Domain%jed, Grid%nk,   &
                             Time%model_time, grid%tracer_axes)
endif  !}

if (do_ocean_age_tracer) then  !{
  call ocean_age_tracer_start(Domain%isd, Domain%ied, Domain%jsd, Domain%jed, T_prog,   &
                              Grid%xt, Grid%yt, Grid%kmt)
endif  !}

#ifdef USE_OCEAN_BGC

if (do_ocmip2_cfc) then  !{
  call ocmip2_cfc_start(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,        &
                        Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                 &
                        T_prog, Time%taup1, Time%model_time,                              &
                        Grid%dat, Grid%tmask, Grid%tracer_axes, Thickness%rho_dzt)
endif  !}

if (do_ocmip2_he) then  !{
  call ocmip2_he_start(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,        &
                       Domain%isd, Domain%ied, Domain%jsd, Domain%jed,                 &
                       T_prog, Time%taup1, Time%model_time,                            &
                       Grid%dat, Grid%tmask, Grid%tracer_axes, Domain%domain2d,        &
                       Thickness%rho_dzt)
endif  !}

if (do_ocean_pert_co2) then  !{
  call ocean_pert_co2_start(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,    &
                            Domain%isd, Domain%ied, Domain%jsd, Domain%jed,             &
                            T_prog, Time%taup1, Time%model_time,                          &
                            grid%dat, grid%tmask,                                       &
                            grid%tracer_axes, Domain%domain2d,                          &
                            Thickness%rho_dzt)
endif  !}

if (do_ocmip2_abiotic) then  !{
  call ocmip2_abiotic_start(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,    &
                            Domain%isd, Domain%ied, Domain%jsd, Domain%jed,             &
                            T_prog, Time%taup1, Time%model_time,                          &
                            Grid%dat, Grid%tmask, Grid%kmt, Grid%xt, Grid%yt, Grid%zt,  &
                            Grid%zw, Grid%dzt,                                          &
                            Grid%tracer_axes, Domain%domain2d,                          &
                            Thickness%rho_dzt)
endif  !}

if (do_ocmip2_biotic) then  !{
  call ocmip2_biotic_start(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,     &
                           Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                           T_prog, Time%taup1, Time%model_time,                           &
                           Grid%dat, Grid%tmask, Grid%kmt, Grid%xt, Grid%yt, Grid%zt,   &
                           Grid%zw, Grid%dzt,                                           &
                           Grid%name, Grid%tracer_axes, Domain%domain2d,                &
                           Thickness%rho_dzt)
endif  !}

if (do_ocean_bgc_restore) then  !{
  call ocean_bgc_restore_start(Domain%isc, Domain%iec, Domain%jsc, Domain%jec,          &
                               Grid%nk, Domain%isd, Domain%ied, Domain%jsd,             &
                               Domain%jed, T_prog, T_diag, Time%taup1, Time%model_time,   &
                               Grid%dat, Grid%tmask, Grid%kmt, Grid%xt, Grid%yt,        &
                               Grid%zt, Grid%zw, Grid%dzt,                              &
                               Grid%name, Grid%tracer_axes, Domain%domain2d,            &
                               Thickness%rho_dzt)
endif  !}

if (do_ocean_po4_pre) then  !{
  call ocean_po4_pre_start(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,     &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,            &
                             T_prog, Time%taup1, Time%model_time,                       &
                             Grid%dat, Grid%tmask, Grid%kmt, Grid%xt, Grid%yt,          &
                             Thickness%rho_dzt)
endif  !}

if (do_ocean_ibgc) then  !{
  call ocean_ibgc_start(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,        &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,            &
                             T_prog, T_diag, Time%taup1, Time%model_time,               &
                             Grid%dat, Grid%tmask, Grid%kmt, Grid%xt, Grid%yt,          &
                             Thickness%depth_zt, Grid%zw, Grid%dzt, Grid%name,&
                             Grid%tracer_axes, Domain%domain2d, Thickness%rho_dzt)

endif  !}

if (do_ocean_fabm) then  !{
  call ocean_fabm_start(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk,     &
                           Domain%isd, Domain%ied, Domain%jsd, Domain%jed,              &
                           T_prog, T_diag, Time%taup1, Time%model_time,                           &
                           Grid%dat, Grid%tmask, Grid%kmt, Grid%xt, Grid%yt, Grid%zt,   &
                           Grid%zw, Grid%dzt,                                           &
                           Grid%name, Grid%tracer_axes, Domain%domain2d,                &
                           Thickness%rho_dzt)
endif  !}

#endif

if (do_transport_matrix) then  !{
  call transport_matrix_start(Time, T_prog, Domain%isd, Domain%ied, Domain%jsd,         &
                              Domain%jed, Grid%nk, Grid%tracer_axes)
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

subroutine ocean_tpm_tracer(Domain, T_prog, T_diag, Grid, Time, Thickness, Dens, dtts, hblt_depth,&
                            sw_pen, opacity, diff_cbt, river)  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

type(ocean_domain_type), intent(in)                             :: Domain
type(ocean_prog_tracer_type), dimension(:), intent(inout)       :: T_prog
type(ocean_diag_tracer_type), dimension(:), intent(inout)       :: T_diag
type(ocean_grid_type), intent(in)                               :: Grid
type(ocean_time_type), intent(in)                               :: Time
type(ocean_thickness_type), intent(in)                          :: Thickness
type(ocean_density_type), intent(in)                            :: Dens
real, intent(in)                                                :: dtts
real, intent(in), dimension(Domain%isd:,Domain%jsd:)            :: hblt_depth
real, intent(in), dimension(Domain%isd:,Domain%jsd:)            :: sw_pen
real, intent(in), dimension(Domain%isd:,Domain%jsd:,:)          :: opacity
real, intent(in), dimension(Domain%isd:,Domain%jsd:,:,:)        :: diff_cbt
real, intent(in), dimension(Domain%isd:,Domain%jsd:)            :: river

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

!
!       set some indices and flags dependent on time
!

call do_time_calc(Time, dtts)

if (do_ocean_age_tracer) then  !{
  call ocean_age_tracer_tracer(T_prog, Time%taup1)
endif  !}

if (do_ocean_residency) then  !{
  call ocean_residency_tracer(Domain%isc, Domain%iec, Domain%jsc, Domain%jec,           &
                              Domain%isd, Domain%ied, Domain%jsd, Domain%jed, Grid%nk,  &
                              T_prog, Grid%tmask, Time%taup1, Time%model_time, dtts)
endif  !}

#ifdef USE_OCEAN_BGC

if (do_ocmip2_abiotic) then  !{
  call ocmip2_abiotic_tracer(Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%nk, &
                             Domain%isd, Domain%ied, Domain%jsd, Domain%jed,          &
                             T_prog, Grid%dat, Grid%tmask, Grid%tcella,               &
                             Time%taum1, dtts, end_of_year)
endif  !}

if (do_ocean_po4_pre) then  !{
  call ocean_po4_pre_tracer(Domain%isc, Domain%iec, Domain%jsc, Domain%jec,          &
                            Domain%isd, Domain%ied, Domain%jsd, Domain%jed, Grid%nk, &
                            T_prog, Time, Thickness, Dens, Grid%zt, hblt_depth)
endif  !}

if (do_ocean_ibgc) then  !{
  call ocean_ibgc_tracer(Domain%isc, Domain%iec, Domain%jsc, Domain%jec,               &
                              Domain%isd, Domain%ied, Domain%jsd, Domain%jed, Grid%nk, &
                              Time, T_prog, Thickness, Dens,                           &
                              Thickness%depth_zt, hblt_depth)

endif  !}

if (do_generic_tracer) then
   call ocean_generic_column_physics(Thickness, hblt_depth, Time, &
        Grid, dtts, Domain%isd,Domain%jsd, T_prog, T_diag,sw_pen,opacity, diff_cbt, Dens )
endif

if (do_ocean_fabm) then  !{
  call ocean_fabm_tracer()

endif  !}

#endif

if (do_transport_matrix) then !{
  call transport_matrix_store_implicit(Time, T_prog, Domain%isd, Domain%ied, Domain%jsd, Domain%jed,    &
       Grid%nk, Domain%isc, Domain%iec, Domain%jsc, Domain%jec, Grid%tmask)
endif  !}


return

end subroutine ocean_tpm_tracer  !}
! </SUBROUTINE> NAME="ocean_tpm_tracer"

end module ocean_tpm_mod  !}

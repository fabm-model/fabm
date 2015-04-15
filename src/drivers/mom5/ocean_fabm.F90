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
#include <fms_platform.h>

!
!
!<CONTACT EMAIL="jorn@bolding-burchard.com"> Jorn Bruggeman
!</CONTACT>
!
!<REVIEWER EMAIL="TODO"> TODO
!</REVIEWER>
!
!<OVERVIEW>
! MOM4 driver for the Framework for Aquatic Biogeochemical Models (FABM)
!</OVERVIEW>
!
!<DESCRIPTION>
!       TODO
!</DESCRIPTION>
!
! <INFO>
! <REFERENCE>
! TODO
! </REFERENCE>
!
! </INFO>
!
! $Id: ocean_fabm.F90 69 2010-09-07 15:00:18Z jornbr $
!

!
!------------------------------------------------------------------
!
!       Module ocean_fabm_mod
!
!       TODO
!
!------------------------------------------------------------------
!

module  ocean_fabm_mod  !{

!
!------------------------------------------------------------------
!
!       Global definitions
!
!------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Modules
!
!----------------------------------------------------------------------
!

use time_manager_mod,         only: time_type
use diag_manager_mod,         only: send_data
use field_manager_mod,        only: fm_field_name_len, fm_path_name_len, fm_string_len, fm_get_index
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value
use fms_mod,                  only: write_data
use mpp_mod,                  only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use mpp_mod,                  only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_manager_mod,         only: get_date
use time_interp_external_mod, only: time_interp_external
use mpp_domains_mod,          only: domain2d, mpp_global_sum, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM
use data_override_mod,only: data_override

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer, otpm_set_diag_tracer
use fm_util_mod, only: fm_util_check_for_bad_fields, fm_util_set_value, fm_util_set_value_string_array
use fm_util_mod, only: fm_util_get_string, fm_util_get_integer_array, fm_util_get_logical, fm_util_get_integer, fm_util_get_real
use fm_util_mod, only: fm_util_get_logical_array, fm_util_get_real_array, fm_util_get_string_array
use fm_util_mod, only: fm_util_start_namelist, fm_util_end_namelist
!use fm_util_mod, only: grid, dtts
!use fm_util_mod, only: taum1, tau, taup1
!use fm_util_mod, only: t_prog, t_diag
use coupler_types_mod,  only: coupler_2d_bc_type ! ind_alpha, ind_csurf
use ocean_types_mod,    only: ocean_prog_tracer_type, ocean_diag_tracer_type, ocean_density_type

!use ocean_sbc_mod,      only: use_waterflux
!use ocean_types_mod,    only: ocean_thickness_type,ocean_density_type

use fabm
use fabm_types
use fabm_config
use fabm_driver

!
!----------------------------------------------------------------------
!
!       force all variables to be "typed"
!
!----------------------------------------------------------------------
!

implicit none

!
!----------------------------------------------------------------------
!
!       Make all routines and variables private by default
!
!----------------------------------------------------------------------
!

private

!
!----------------------------------------------------------------------
!
!       Public routines
!
!----------------------------------------------------------------------
!

public  :: ocean_fabm_bbc
public  :: ocean_fabm_end
public  :: ocean_fabm_init
public  :: ocean_fabm_sbc
public  :: ocean_fabm_source
public  :: ocean_fabm_start
public  :: ocean_fabm_tracer
public  :: ocean_fabm_avg_sfc
public  :: ocean_fabm_sum_sfc
public  :: ocean_fabm_zero_sfc
public  :: ocean_fabm_flux_init
public  :: ocean_fabm_sfc_end
public  :: ocean_fabm_init_sfc

!
!----------------------------------------------------------------------
!
!       Private routines
!
!----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Private parameters
!
!----------------------------------------------------------------------
!

character(len=fm_field_name_len), parameter     :: package_name = 'ocean_fabm'
character(len=48), parameter                    :: mod_name = 'ocean_fabm_mod'
character(len=fm_string_len), parameter         :: default_restart_file = 'ocean_fabm.res.nc'
character(len=fm_string_len), parameter         :: default_ice_restart_file = 'ice_ocean_fabm.res.nc'
character(len=fm_string_len), parameter         :: default_ocean_restart_file = 'ocean_fabm_airsea_flux.res.nc'

!
!----------------------------------------------------------------------
!
!       Private types
!
!----------------------------------------------------------------------
!


type biotic_type  !{

  type (type_model),pointer :: model
  type (type_bulk_variable_id)                     :: id_temp,id_salt,id_pres,id_par,id_dens

  character(len=fm_field_name_len)                 :: name
  logical                                          :: do_virtual_flux = .false.
  integer,_ALLOCATABLE,dimension(:)                :: inds,inds_clip,inds_diag,inds_diag_hz,inds_cons,inds_cons_tot,inds_cons_ave
  double precision,_ALLOCATABLE,dimension(:,:,:,:) :: work_cons _NULL
  double precision,_ALLOCATABLE,dimension(:,:,:)   :: w _NULL,sf_fluxes _NULL
  double precision,_ALLOCATABLE,dimension(:,:)     :: adv _NULL,work_dy _NULL

  real,         _ALLOCATABLE               :: bottom_state(:,:,:) _NULL
  real,         _ALLOCATABLE               :: surface_state(:,:,:) _NULL

  ! Arrays to hold information on externally provided fields
  character*128,_ALLOCATABLE, dimension(:) :: ext_2d_variables _NULL,ext_3d_variables _NULL
  real,         _ALLOCATABLE               :: ext_2d_data(:,:,:) _NULL, ext_3d_data(:,:,:,:) _NULL

  integer :: ind_diag_clip
end type biotic_type  !}

   type,extends(type_base_driver) :: type_mom_driver
   contains
      procedure :: mom_driver_fatal_error
      procedure :: mom_driver_log_message
   end type

double precision,_ALLOCATABLE,dimension(:,:,:) :: clipped _NULL

!
!----------------------------------------------------------------------
!
!       Public variables
!
!----------------------------------------------------------------------
!

logical, public :: do_ocean_fabm

!
!----------------------------------------------------------------------
!
!       Private variables
!
!----------------------------------------------------------------------
!

integer                                 :: package_index
logical                                 :: module_initialized = .false.
integer                                 :: indsal
integer                                 :: indtemp

integer                                 :: id_clock_fabm_source,id_clock_fabm_vertmov,id_clock_fabm_sbc

!
!----------------------------------------------------------------------
!
!       Input parameters:
!
!----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Calculated parameters (with possible initial input values):
!
!----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Global atmosphere model
!
!----------------------------------------------------------------------
!
integer                                             :: instances
type(biotic_type), allocatable, dimension(:),target :: biotic
integer                                             :: index_irr,index_chl
logical                                             :: inplace_repair,zero_river_concentration
logical                                             :: disable_sources,disable_vertical_movement

!
!-----------------------------------------------------------------------
!
!       Subroutine and function definitions
!
!-----------------------------------------------------------------------
!

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocean_fabm_init  !{

use fms_mod, only : open_namelist_file,close_file,file_exist
use mpp_io_mod, only: mpp_open,MPP_RDONLY
use diag_manager_mod, only: register_diag_field

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'


!
!-----------------------------------------------------------------------
!       arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=fm_string_len)                            :: string
character(len=fm_field_name_len+3)                      :: long_suffix
character(len=256)                                      :: caller_str
character(len=fm_string_len), pointer, dimension(:)     :: good_list
integer                                                 :: n,extcount

!integer                                                 :: n
!character(len=fm_field_name_len+1)                      :: suffix
!character(len=fm_field_name_len+3)                      :: long_suffix
!character(len=fm_field_name_len)                        :: name
integer                                                 :: nmlunit
integer                                                 :: i
character(len=256)                                      :: namelist_file

!
!       Initialize the FABM package
!

package_index = otpm_set_tracer_package(package_name,            &
     restart_file = default_restart_file,     &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')

!
!       Check whether to use this package
!

path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances < 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
endif  !}

!
!       Check some things
!

write (stdout(),*)
if (instances == 0) then  !{
  write (stdout(),*) trim(note_header), ' No instances'
  do_ocean_fabm = .false.
else  !}{
  if (instances == 1) then  !{
    write (stdout(),*) trim(note_header), ' ', instances, ' instance'
  else  !}{
    write (stdout(),*) trim(note_header), ' ', instances, ' instances'
  endif  !}
  do_ocean_fabm = .true.
endif  !}

module_initialized = .true.

!
!       Return if we don't want to use this package,
!       after changing the list back
!

if (.not. do_ocean_fabm) then  !{
  return
endif  !}

! Connect custom logging/error reporting object to FABM.
allocate(type_mom_driver::driver)

!
!       allocate storage for biotic array
!

allocate ( biotic(instances) )

!
!       loop over the names, saving them into the biotic array
!

do n = 1, instances  !{

  if (fm_get_value(path_to_names, name, index = n)) then  !{
    biotic(n)%name = name
  else  !}{
    write (name,*) n
    call mpp_error(FATAL, trim(error_header) //        &
         'Bad field name for index ' // trim(name))
  endif  !}

enddo  !}

!
!-----------------------------------------------------------------------
!       Process the namelists
!-----------------------------------------------------------------------
!

!
!       Add the package name to the list of good namelists, to be used
!       later for a consistency check
!

if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) <= 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif  !}

!
!-----------------------------------------------------------------------
!       Set up the *global* namelist (shared across all instances)
!-----------------------------------------------------------------------
!

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

call fm_util_start_namelist(package_name, '*global*', caller = caller_str, no_overwrite = .true., &
     check = .true.)

call fm_util_set_value('wind_file', 'INPUT/scalar_wind_ongrid.nc')
call fm_util_set_value('wind_name', 'scalar_wind')
call fm_util_set_value('inplace_repair', .false.)
call fm_util_set_value('zero_river_concentration', .false.)
call fm_util_set_value('disable_sources', .false.)
call fm_util_set_value('disable_vertical_movement', .false.)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)

!
!-----------------------------------------------------------------------
!       Set up the instance namelists
!-----------------------------------------------------------------------
!

!do n = 1, instances  !{
!
!!
!!       create the instance namelist
!!
!
!  call fm_util_start_namelist(package_name, biotic(n)%name, caller = caller_str, no_overwrite = .true., &
!       check = .true.)
!
!  call fm_util_set_value('namelist_file', 'fabm.nml')
!
!  call fm_util_end_namelist(package_name, biotic(n)%name, check = .true., caller = caller_str)
!
!enddo  !} n


do n = 1, instances  !{

  call fm_util_start_namelist(package_name, biotic(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)
  call fm_util_set_value('namelist_file', 'fabm.nml')

!  call fm_util_start_namelist(package_name, biotic(n)%name, caller = caller_str)
  namelist_file = fm_util_get_string ('namelist_file', caller = caller_str, scalar = .true.)

  name = biotic(n)%name
  if (name(1:1) == '_') then  !{
    suffix = ' '
    long_suffix = ' '
  else  !}{
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif  !}

  ! Create the FABM model tree
  if (file_exist('fabm.yaml')) then
      call mpp_open(nmlunit, 'fabm.yaml', action=MPP_RDONLY )
      allocate(biotic(n)%model)
      call fabm_create_model_from_yaml_file(biotic(n)%model,unit=nmlunit)
  else
      nmlunit = open_namelist_file(namelist_file)
      biotic(n)%model => fabm_create_model_from_file(nmlunit)
  end if
  call close_file (nmlunit)

  ! Register state variables
  allocate(biotic(n)%inds(size(biotic(n)%model%info%state_variables)))
  do i=1,size(biotic(n)%model%info%state_variables)
     biotic(n)%inds(i) = otpm_set_prog_tracer(                                               &
          trim(biotic(n)%model%info%state_variables(i)%name) // suffix,                      &
          package_name,                                                                      &
          longname = trim(biotic(n)%model%info%state_variables(i)%long_name) // trim(long_suffix),  &
          units = trim(biotic(n)%model%info%state_variables(i)%units),                       &
          flux_units = trim(biotic(n)%model%info%state_variables(i)%units)//'/s',            &
          caller = trim(mod_name)//'('//trim(sub_name)//')',                                 &
          min_tracer_limit = biotic(n)%model%info%state_variables(i)%minimum,                &
          max_tracer_limit = biotic(n)%model%info%state_variables(i)%maximum,                &
          !min_tracer = biotic(n)%model%info%state_variables(i)%minimum,                      &
          !max_tracer = biotic(n)%model%info%state_variables(i)%maximum,                      &
          const_init_tracer = .true.,                                                        &
          const_init_value = biotic(n)%model%info%state_variables(i)%initial_value)
  end do

!  call fm_util_end_namelist(package_name, biotic(n)%name, caller = caller_str)
  call fm_util_end_namelist(package_name, biotic(n)%name, check = .true., caller = caller_str)

  ! Obtain ids of required external variables
  biotic(n)%id_temp   = fabm_get_bulk_variable_id(biotic(n)%model,standard_variables%temperature)
  biotic(n)%id_salt   = fabm_get_bulk_variable_id(biotic(n)%model,standard_variables%practical_salinity)
  biotic(n)%id_pres   = fabm_get_bulk_variable_id(biotic(n)%model,standard_variables%pressure)
  biotic(n)%id_dens   = fabm_get_bulk_variable_id(biotic(n)%model,standard_variables%density)
  biotic(n)%id_par    = fabm_get_bulk_variable_id(biotic(n)%model,standard_variables%downwelling_photosynthetic_radiative_flux)

enddo  !} n

!
!       Check for any errors in the number of fields in the namelists for this package
!

good_list => fm_util_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',   &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
if (associated(good_list)) then  !{
  call fm_util_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,       &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
endif  !}

index_chl = otpm_set_diag_tracer('chl',                                             &
 caller=trim(mod_name)//'('//trim(sub_name)//')',                                      &
 longname='chlorophyll', units='mg/m^3',                            &
 conversion=1.0, offset=0.0, min_tracer=0.0, max_tracer=1.e20, &
 min_range=-10.0, max_range=100.0, const_init_tracer=.true.,const_init_value=0.0,      &
 restart_file='ocean_chl.res.nc')

index_irr = otpm_set_diag_tracer('irr',                                             &
 caller=trim(mod_name)//'('//trim(sub_name)//')',                                      &
 longname='irradiance', units='W/m^2',                        &
 conversion=1.0, offset=0.0, min_tracer=0.0, max_tracer=1.e20, &
 min_range=-10.0, max_range=100.0, const_init_tracer=.true.,const_init_value=0.0,      &
 restart_file='ocean_irr.res.nc')

if (index_irr==-1 .or. index_chl==-1) then
    call mpp_error(FATAL,trim(error_header) // 'index of chlorophyll and/or irradiance equals -1.')
end if

id_clock_fabm_source = mpp_clock_id('(Ocean FABM: source) ',grain=CLOCK_ROUTINE)
id_clock_fabm_vertmov = mpp_clock_id('(Ocean FABM: vertical movement) ',grain=CLOCK_ROUTINE)
id_clock_fabm_sbc = mpp_clock_id('(Ocean FABM: sbc) ',grain=CLOCK_ROUTINE)

return

99 call mpp_error(FATAL,trim(error_header) // ' Unable to open FABM namelist file "' // trim(namelist_file) // '"')

end subroutine ocean_fabm_init  !}
! </SUBROUTINE> NAME="ocean_fabm_init"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!

subroutine ocean_fabm_start(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,      &
     T_prog, T_diag, taup1, model_time, grid_dat, grid_tmask, grid_kmt,                 &
     grid_xt, grid_yt, grid_zt, grid_zw, grid_dzt, grid_name, grid_tracer_axes, &
     mpp_domain2d, rho_dzt, dzt, swflx, swflx_vis, current_wave_stress)  !{

!
!-----------------------------------------------------------------------
!       modules (have to come first)
!-----------------------------------------------------------------------
!

use time_interp_external_mod, only: init_external_field
use diag_manager_mod, only        : register_diag_field, diag_axis_init
use field_manager_mod, only       : fm_get_index

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
type(ocean_diag_tracer_type), intent(in), dimension(:)  :: T_diag
integer, intent(in)                                     :: taup1
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
integer, dimension(isd:,jsd:), intent(in)               :: grid_kmt
real, dimension(isd:,jsd:), intent(in)                  :: grid_xt
real, dimension(isd:,jsd:), intent(in)                  :: grid_yt
real, dimension(:), intent(in)                          :: grid_zt
real, dimension(:), intent(in)                          :: grid_zw
real, dimension(:), intent(in)                          :: grid_dzt
character(len=*), intent(in)                            :: grid_name
integer, dimension(3), intent(in)                       :: grid_tracer_axes
type(domain2d), intent(in)                              :: mpp_domain2d
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt
real, dimension(isd:,jsd:,:), target                    :: dzt
real, dimension(isd:,jsd:), target                      :: swflx, swflx_vis, current_wave_stress

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_start'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!
character(len=256)                                      :: caller_str
character(len=fm_field_name_len+1)                      :: str
integer                                                 :: i,n=1,index_wind

!
! =====================================================================
!       begin of executable code
! =====================================================================
!
!
!-----------------------------------------------------------------------
!       give info
!-----------------------------------------------------------------------
!

write(stdout(),*)
write(stdout(),*) trim(note_header),                     &
                  'Starting ', trim(package_name), ' module'

indtemp = fm_get_index('/ocean_mod/prog_tracers/temp')
if (indtemp <= 0) then  !{
  call mpp_error(FATAL,trim(error_header) // ' Could not get the temperature index')
endif  !}

indsal = fm_get_index('/ocean_mod/prog_tracers/salt')
if (indsal <= 0) then  !{
  call mpp_error(FATAL,trim(error_header) // ' Could not get the salinity index')
endif  !}
!
!-----------------------------------------------------------------------
!       save the *global* namelist values
!-----------------------------------------------------------------------
!

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

call fm_util_start_namelist(package_name, '*global*', caller = caller_str)

inplace_repair            = fm_util_get_logical('inplace_repair',           scalar = .true.)
zero_river_concentration  = fm_util_get_logical('zero_river_concentration', scalar = .true.)
disable_sources           = fm_util_get_logical('disable_sources',          scalar = .true.)
disable_vertical_movement = fm_util_get_logical('disable_vertical_movement',scalar = .true.)

call fm_util_end_namelist(package_name, '*global*', caller = caller_str)

!
!-----------------------------------------------------------------------
!       read in the namelists for each instance
!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!       Set up diagnostic fields
!-----------------------------------------------------------------------
do n = 1, instances  !{

  if (instances == 1) then  !{
    str = ' '
  else  !}{
    str = '_' // biotic(n)%name
  endif  !}

  ! Register diagnostic variables defined on full model domain.
  allocate(biotic(n)%inds_diag(size(biotic(n)%model%info%diagnostic_variables)))
  do i=1,size(biotic(n)%model%info%diagnostic_variables)
     biotic(n)%inds_diag(i) = register_diag_field('ocean_model',      &
          trim(biotic(n)%model%info%diagnostic_variables(i)%name)//str, grid_tracer_axes(1:3),                       &
          model_time, trim(biotic(n)%model%info%diagnostic_variables(i)%long_name), &
          trim(biotic(n)%model%info%diagnostic_variables(i)%units),            &
          missing_value = biotic(n)%model%info%diagnostic_variables(i)%missing_value)
     biotic(n)%model%info%diagnostic_variables(i)%save = biotic(n)%inds_diag(i)/=-1
  end do

  ! Register diagnostic variables defined on horizontal slice of domain.
  allocate(biotic(n)%inds_diag_hz(size(biotic(n)%model%info%diagnostic_variables_hz)))
  do i=1,size(biotic(n)%model%info%diagnostic_variables_hz)
     biotic(n)%inds_diag_hz(i) = register_diag_field('ocean_model',      &
          trim(biotic(n)%model%info%diagnostic_variables_hz(i)%name)//str, grid_tracer_axes(1:2),                       &
          model_time, trim(biotic(n)%model%info%diagnostic_variables_hz(i)%long_name), &
          trim(biotic(n)%model%info%diagnostic_variables_hz(i)%units),            &
          missing_value = biotic(n)%model%info%diagnostic_variables_hz(i)%missing_value)
     biotic(n)%model%info%diagnostic_variables_hz(i)%save = biotic(n)%inds_diag_hz(i)/=-1
  end do

  ! Register clipping diagnostic (increase/time) for pelagic state variables.
  allocate(biotic(n)%inds_clip(size(biotic(n)%model%info%state_variables)))
  do i=1,size(biotic(n)%model%info%state_variables)
     biotic(n)%inds_clip(i) = register_diag_field('ocean_model',      &
          trim(biotic(n)%model%info%state_variables(i)%name)//'_clip'//str, grid_tracer_axes(1:3),                       &
          model_time, trim(biotic(n)%model%info%state_variables(i)%long_name)//' clipping increase', &
          trim(biotic(n)%model%info%state_variables(i)%units),            &
          missing_value = -1.0e+10)
  end do

  ! Register local values anf domain-wide integrals of conserved quantities.
  allocate(biotic(n)%inds_cons    (size(biotic(n)%model%info%conserved_quantities)))
  allocate(biotic(n)%inds_cons_tot(size(biotic(n)%model%info%conserved_quantities)))
  allocate(biotic(n)%inds_cons_ave(size(biotic(n)%model%info%conserved_quantities)))
  do i=1,size(biotic(n)%model%info%conserved_quantities)
     biotic(n)%inds_cons(i) = register_diag_field('ocean_model',                            &
          trim(biotic(n)%model%info%conserved_quantities(i)%name)//str, grid_tracer_axes(1:3),                 &
          model_time, trim(biotic(n)%model%info%conserved_quantities(i)%long_name), &
          trim(biotic(n)%model%info%conserved_quantities(i)%units),                     &
          missing_value = -1.0e+10)
     biotic(n)%inds_cons_tot(i) = register_diag_field('ocean_model',                            &
          trim(biotic(n)%model%info%conserved_quantities(i)%name)//'_total'//str,                 &
          model_time, 'mass integrated '//trim(biotic(n)%model%info%conserved_quantities(i)%long_name), &
          trim(biotic(n)%model%info%conserved_quantities(i)%units)//'*m3/1e21',                     &
          missing_value = -1.0e+10)
     biotic(n)%inds_cons_ave(i) = register_diag_field('ocean_model',                            &
          trim(biotic(n)%model%info%conserved_quantities(i)%name)//'_global_ave'//str,                 &
          model_time, 'global mass weighted mean '//trim(biotic(n)%model%info%conserved_quantities(i)%long_name)//' in liquid seawater', &
          trim(biotic(n)%model%info%conserved_quantities(i)%units),                     &
          missing_value = -1.0e+10)
  end do

  biotic(n)%ind_diag_clip = register_diag_field('ocean_model',      &
    'fabm_clip', grid_tracer_axes(1:3), model_time, 'fabm_clip', &
    '-', missing_value = -1.0e+10)

  call fabm_set_domain(biotic(n)%model,iec-isc+1,jec-jsc+1,nk)
  call fabm_set_mask(biotic(n)%model,grid_tmask(isc:iec,jsc:jec,:))

enddo  !} n

!
!-----------------------------------------------------------------------
!     give info
!-----------------------------------------------------------------------
!

write(stdout(),*)
write(stdout(),*) trim(note_header), 'Tracer runs initialized'
write(stdout(),*)

allocate(clipped(isc:iec,jsc:jec,nk))

do n=1,instances

   allocate(biotic(n)%work_dy(isc:iec,     size(biotic(n)%model%info%state_variables)))
   allocate(biotic(n)%w      (isc:iec,nk+1,size(biotic(n)%model%info%state_variables)))
   allocate(biotic(n)%adv    (        nk+1,size(biotic(n)%model%info%state_variables)))
   if (any(biotic(n)%inds_cons>0).or.any(biotic(n)%inds_cons_tot>0).or.any(biotic(n)%inds_cons_ave>0)) then
      allocate(biotic(n)%work_cons(isd:ied,jsd:jed,nk,size(biotic(n)%model%info%conserved_quantities)))
      biotic(n)%work_cons = 0
   end if
  allocate(biotic(n)%bottom_state(isc:iec,jsc:jec,size(biotic(n)%model%info%bottom_state_variables)))
  do i=1,size(biotic(n)%model%info%bottom_state_variables)
      call fabm_link_bottom_state_data(biotic(n)%model,i,biotic(n)%bottom_state(:,:,i))
  end do

   !if (biotic(n)%id_pres  /=-1) call fabm_link_bulk_data   (biotic(n)%model,biotic(n)%id_pres,  Dens%pressure_at_depth (isc:iec,jsc:jec,:))
   if (fabm_variable_needs_values(biotic(n)%model,biotic(n)%id_par)) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_par,  t_diag(index_irr)%field(isc:iec,jsc:jec,:))

   call fabm_link_horizontal_data(biotic(n)%model,standard_variables%surface_downwelling_shortwave_flux,swflx(isc:iec,jsc:jec))
   call fabm_link_horizontal_data(biotic(n)%model,standard_variables%surface_downwelling_photosynthetic_radiative_flux,swflx_vis(isc:iec,jsc:jec))
   call fabm_link_horizontal_data(biotic(n)%model,standard_variables%bottom_stress,current_wave_stress(isc:iec,jsc:jec))
   call fabm_link_bulk_data(biotic(n)%model,standard_variables%cell_thickness,dzt(isc:iec,jsc:jec,:))
end do

return

end subroutine  ocean_fabm_start  !}
! </SUBROUTINE> NAME="ocean_fabm_start"

subroutine get_external_fields(isc, iec, jsc, jec, nk, model_time)

integer,intent(in)                                      :: isc, iec, jsc, jec, nk
type(time_type), intent(in)                             :: model_time
real, allocatable                                       :: dummy(:,:,:),dummy_hz(:,:)
logical                                                 :: used
integer                                                 :: extcount,n,i

do n = 1, instances  !{
   if (.not. _ALLOCATED(biotic(n)%ext_3d_variables)) then
      ! This is the first time that this subroutine is called. Therefore, probe for
      ! each FABM variable whether external fields have been provided. If so, register
      ! the field name and increment the number of external fields. The total count
      ! is then used to allocate arrays for the external data, which are then coupled to
      ! FABM.
      !
      ! It would have been nice to do this during initialization (_init or _start),
      ! instead of on demand, but during initialization the data override module has not
      ! been told the ocean domain yet. Thus, data_override cannot be called.

      ! Process 3D variables.
      extcount = 0
      allocate(biotic(n)%ext_3d_variables(1:size(biotic(n)%model%info%dependencies)))
      biotic(n)%ext_3d_variables = ''
      allocate(dummy(isc:iec,jsc:jec,1:nk))  ! 3D dummy array used to receive external values, needed as long as the real target array is not allocated yet
      do i=1,size(biotic(n)%model%info%dependencies)
         call data_override ('OCN', 'fabm_'//trim(biotic(n)%model%info%dependencies(i)), dummy, model_time, override=used )
         if (used) then
            extcount = extcount + 1
            biotic(n)%ext_3d_variables(extcount) = trim(biotic(n)%model%info%dependencies(i))
         end if
      end do
      deallocate(dummy)
      allocate(biotic(n)%ext_3d_data(isc:iec,jsc:jec,1:nk,1:extcount))
      do i=1,size(biotic(n)%ext_3d_variables)
         if (biotic(n)%ext_3d_variables(i)=='') exit
         call fabm_link_bulk_data(biotic(n)%model,biotic(n)%ext_3d_variables(i),biotic(n)%ext_3d_data(isc:iec,jsc:jec,1:nk,i))
      end do

      ! Process 2D variables.
      extcount = 0
      allocate(biotic(n)%ext_2d_variables(1:size(biotic(n)%model%info%dependencies_hz)))
      biotic(n)%ext_2d_variables = ''
      allocate(dummy_hz(isc:iec,jsc:jec))  ! 2D dummy array used to receive external values, needed as long as the real target array is not allocated yet
      do i=1,size(biotic(n)%model%info%dependencies_hz)
         call data_override ('OCN', 'fabm_'//trim(biotic(n)%model%info%dependencies_hz(i)), dummy_hz, model_time, override=used )
         if (used) then
            extcount = extcount + 1
            biotic(n)%ext_2d_variables(extcount) = trim(biotic(n)%model%info%dependencies_hz(i))
         end if
      end do
      deallocate(dummy_hz)
      allocate(biotic(n)%ext_2d_data(isc:iec,jsc:jec,1:extcount))
      do i=1,size(biotic(n)%ext_2d_variables)
         if (biotic(n)%ext_2d_variables(i)=='') exit
         call fabm_link_horizontal_data(biotic(n)%model,biotic(n)%ext_2d_variables(i),biotic(n)%ext_2d_data(isc:iec,jsc:jec,i))
      end do
   end if

   ! Update values of external variables.
   do i=1,ubound(biotic(n)%ext_3d_data,4)
      call data_override('OCN','fabm_'//trim(biotic(n)%ext_3d_variables(i)), biotic(n)%ext_3d_data(:,:,:,i), model_time)
   end do
   do i=1,ubound(biotic(n)%ext_2d_data,3)
      call data_override('OCN','fabm_'//trim(biotic(n)%ext_2d_variables(i)), biotic(n)%ext_2d_data(:,:,i), model_time)
   end do
end do

end subroutine get_external_fields

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_source">
!
! <DESCRIPTION>
!     compute the source terms for the BIOTICs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!

subroutine ocean_fabm_source(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,     &
     T_prog, T_diag, taum1, model_time, grid_dat, grid_zw, grid_ht, grid_tmask, grid_kmt, rho_dzt, dzt, Dens, dtts, mpp_domain2d)  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: nk
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
type(ocean_diag_tracer_type), intent(inout), dimension(:)       :: T_diag
integer, intent(in)                                             :: taum1
type(time_type), intent(in)                                     :: model_time
real, dimension(isd:,jsd:), intent(in)                          :: grid_dat
real, dimension(nk), intent(in)                                 :: grid_zw
real, dimension(isd:,jsd:), intent(in)                          :: grid_ht
real, dimension(isd:,jsd:,:), intent(in)                        :: grid_tmask
integer, dimension(isd:,jsd:), intent(in)                       :: grid_kmt
real, dimension(isd:,jsd:,:,:), intent(in)                      :: rho_dzt
real, dimension(isd:ied,jsd:jed,nk), intent(in)                 :: dzt
type(ocean_density_type), intent(in)                            :: Dens
real, intent(in)                                                :: dtts
type(domain2d), intent(in)                                      :: mpp_domain2d

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_source'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i
integer :: j
integer :: k
integer :: n
integer :: ivar
logical :: used,valid
real, dimension(isd:ied,jsd:jed) :: tracer_k
real :: total_tracer,total_volume

!
! =====================================================================
!     begin executable code
! =====================================================================
!

call mpp_clock_begin(id_clock_fabm_source)

! Set a reasonable chlorophyll concentration
T_diag(index_chl)%field(isc:iec,jsc:jec,:) = 0.08

! Update values of external variables.
call get_external_fields(isc, iec, jsc, jec, nk, model_time)

do n = 1, instances  !{

  ! Set pointers to environmental variables.
  if (fabm_is_variable_used(biotic(n)%id_temp)) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_temp,t_prog(indtemp)%field(isc:iec,jsc:jec,:,taum1))
  if (fabm_is_variable_used(biotic(n)%id_salt)) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_salt,t_prog(indsal )%field(isc:iec,jsc:jec,:,taum1))
  if (fabm_is_variable_used(biotic(n)%id_dens)) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_dens,Dens%rho             (isc:iec,jsc:jec,:,taum1))
  if (fabm_is_variable_used(biotic(n)%id_pres)) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_pres,Dens%pressure_at_depth(isc:iec,jsc:jec,:))
  call fabm_link_horizontal_data(biotic(n)%model,standard_variables%bottom_depth_below_geoid,grid_ht(isc:iec,jsc:jec))

  ! Link to current biogeochemical state variable values, as maintained by MOM4.
  do ivar=1,size(biotic(n)%model%info%state_variables)
     t_prog(biotic(n)%inds(ivar))%wrk1(isc:iec,jsc:jec,:) = t_prog(biotic(n)%inds(ivar))%field(isc:iec,jsc:jec,:,taum1)
     call fabm_link_bulk_state_data(biotic(n)%model,ivar,t_prog(biotic(n)%inds(ivar))%wrk1(isc:iec,jsc:jec,:))
  end do

  call fabm_check_ready(biotic(n)%model)

  ! Repair bio state at the start of the time step
  clipped = 0
  do k = 1, nk  !{
    do j = jsc, jec  !{
      call fabm_check_state(biotic(n)%model,1,iec-isc+1,j-jsc+1,k,.true.,valid)
    end do
  end do

  ! Send per-grid-point clipping status (0 = no clipping, 1 = clipped) to diagnostic manager
  if (biotic(n)%ind_diag_clip>0) &
     used = send_data(biotic(n)%ind_diag_clip,clipped(isc:iec,jsc:jec,:),model_time,rmask=grid_tmask(isc:iec,jsc:jec,:))

  ! Send per-grid-point, per-variable clipping-induced change to diagnostic manager
  do ivar=1,size(biotic(n)%model%info%state_variables)
     if (biotic(n)%inds_clip(ivar)>0) then
        clipped(isc:iec,jsc:jec,:) = t_prog(biotic(n)%inds(ivar))%wrk1(isc:iec,jsc:jec,:) - t_prog(biotic(n)%inds(ivar))%field(isc:iec,jsc:jec,:,taum1)
        used = send_data(biotic(n)%inds_clip(ivar),clipped(isc:iec,jsc:jec,:),model_time,rmask=grid_tmask(isc:iec,jsc:jec,:))
     end if
  end do

  if (any(biotic(n)%inds_cons>0).or.any(biotic(n)%inds_cons_tot>0).or.any(biotic(n)%inds_cons_ave>0)) then
     ! Get 3D fields of conserved quantities from FABM
     do k = 1, nk  !{
       do j = jsc, jec  !{
         call fabm_get_conserved_quantities(biotic(n)%model,1,iec-isc+1,j-jsc+1,k,biotic(n)%work_cons(isc:iec,j,k,:))
       end do
     end do

     ! If we need the global average for any variable, calculate the global integral of mass here.
     if (any(biotic(n)%inds_cons_ave>0)) then
        total_volume = 0.0
        do k=1,nk
           tracer_k(:,:) =  grid_tmask(:,:,k)*grid_dat(:,:)*rho_dzt(:,:,k,taum1) !*dzt(:,:,k)
           !if(have_obc) tracer_k(:,:) = tracer_k(:,:)*Grd%obc_tmask(:,:)
           total_volume  = total_volume + mpp_global_sum(mpp_domain2d,tracer_k(:,:), NON_BITWISE_EXACT_SUM)
        end do
     end if

     do ivar=1,size(biotic(n)%model%info%conserved_quantities)
        ! Send full 3D conserved quantity field to diagnostic manager.
        if (biotic(n)%inds_cons(ivar)>0) then
           used = send_data(biotic(n)%inds_cons(ivar),biotic(n)%work_cons(isc:iec,jsc:jec,:,ivar),model_time,rmask=grid_tmask(isc:iec,jsc:jec,:))
        end if

        if (biotic(n)%inds_cons_tot(ivar)>0.or.biotic(n)%inds_cons_ave(ivar)>0) then
           ! Integrate tracer across full domain.
           total_tracer = 0.0
           do k=1,nk
              tracer_k(:,:) =  grid_tmask(:,:,k)*grid_dat(:,:)*biotic(n)%work_cons(:,:,k,ivar)*rho_dzt(:,:,k,taum1) !*dzt(:,:,k)
              !if(have_obc) tracer_k(:,:) = tracer_k(:,:)*Grd%obc_tmask(:,:)
              total_tracer  = total_tracer + mpp_global_sum(mpp_domain2d,tracer_k(:,:), NON_BITWISE_EXACT_SUM)
           end do

           ! Send global integral and mean to diagnostic manager.
           if (biotic(n)%inds_cons_tot(ivar)>0) used = send_data(biotic(n)%inds_cons_tot(ivar),total_tracer*1e-21,model_time)
           if (biotic(n)%inds_cons_ave(ivar)>0) used = send_data(biotic(n)%inds_cons_ave(ivar),total_tracer/total_volume,model_time)
        end if
     end do
  end if

end do

!
!       Loop over multiple instances
!
do n = 1, instances  !{

  do k = 1, nk  !{
    do j = jsc, jec  !{

      ! Initialize derivatives to zero, because the bio model will increment/decrement values rather than set them.
      biotic(n)%work_dy = 0.d0

      call fabm_do(biotic(n)%model,1,iec-isc+1,j-jsc+1,k,biotic(n)%work_dy(isc:iec,:))
      if (any(isnan(biotic(n)%work_dy(isc:iec,:)))) then
         call mpp_error(FATAL,trim(error_header) // ' NaN in FABM sink/source terms.')
      end if

      if (.not.disable_sources) then
         ! Update tendencies with current sink and source terms.
         do ivar=1,size(biotic(n)%model%info%state_variables)
            if (inplace_repair) then
               ! Add clipping difference divided by time step as tracer source term
               biotic(n)%work_dy(isc:iec,ivar) = biotic(n)%work_dy(isc:iec,ivar) + &
                                     (t_prog(biotic(n)%inds(ivar))%wrk1(isc:iec,j,k) - t_prog(biotic(n)%inds(ivar))%field(isc:iec,j,k,taum1))/dtts
            endif

            t_prog(biotic(n)%inds(ivar))%th_tendency(isc:iec,j,k) = t_prog(biotic(n)%inds(ivar))%th_tendency(isc:iec,j,k) &
               + biotic(n)%work_dy(isc:iec,ivar)*rho_dzt(isc:iec,j,k,taum1)*grid_tmask(isc:iec,j,k)
         end do
      end if

    end do  !} j
  end do  !} k

  ! Save diagnostic variables defined on full domain.
  do ivar=1,size(biotic(n)%model%info%diagnostic_variables)
      if (biotic(n)%inds_diag(ivar) > 0) then
        used = send_data(biotic(n)%inds_diag(ivar),               &
             fabm_get_bulk_diagnostic_data(biotic(n)%model,ivar), &
             model_time, rmask = grid_tmask(isc:iec,jsc:jec,:))
      endif
  end do

  ! Save diagnostic variables defined on horizontal domain only.
  do ivar=1,size(biotic(n)%model%info%diagnostic_variables_hz)
      if (biotic(n)%inds_diag_hz(ivar) > 0) then
        used = send_data(biotic(n)%inds_diag_hz(ivar),                  &
             fabm_get_horizontal_diagnostic_data(biotic(n)%model,ivar), &
             model_time, rmask = grid_tmask(isc:iec,jsc:jec,1))
      endif
  end do
end do

call mpp_clock_end(id_clock_fabm_source)

call mpp_clock_begin(id_clock_fabm_vertmov)

do n = 1, instances  !{

  ! Vertical movement is applied with a first-order upwind scheme.
  do j = jsc, jec
    ! For every i: get sinking speed in cell centers, over entire column, for all state variables.
    do k=1,nk
      call fabm_get_vertical_movement(biotic(n)%model,1,iec-isc+1,j-jsc+1,k,biotic(n)%w(:,k,:))
    end do

    do i = isc, iec
     if (grid_tmask(i,j,1)/=1.) cycle

     ! Interpolate to sinking speed at interfaces
     do ivar=1,size(biotic(n)%model%info%state_variables)
       biotic(n)%w(i,2:grid_kmt(i,j),ivar) = (biotic(n)%w(i,2:grid_kmt(i,j),ivar) + biotic(n)%w(i,1:grid_kmt(i,j)-1,ivar))*0.5d0
     end do
     biotic(n)%w(i,1,              :) = 0.0d0   ! Surface boundary condition
     biotic(n)%w(i,grid_kmt(i,j)+1,:) = 0.0d0   ! Bottom boundary condition

     biotic(n)%adv = 0.d0

     do ivar=1,size(biotic(n)%model%info%state_variables)

        ! Get tracer flux at all interfaces.
        do k=2,grid_kmt(i,j)
          if (biotic(n)%w(i,k,ivar)>0.) then
            ! floating
            biotic(n)%adv(k,ivar) = biotic(n)%w(i,k,ivar)*t_prog(biotic(n)%inds(ivar))%field(i,j,k,taum1)*Dens%rho(i,j,k,taum1)
          elseif (biotic(n)%w(i,k,ivar)<0.) then
            ! sinking
            biotic(n)%adv(k,ivar) = biotic(n)%w(i,k,ivar)*t_prog(biotic(n)%inds(ivar))%field(i,j,k-1,taum1)*Dens%rho(i,j,k-1,taum1)
          end if
        end do

        if (.not.disable_vertical_movement) then
           ! Add transport to tracer tendencies (conservative formulation)
           ! Note: normally the transport terms should be divided by the layer thickness.
           ! However, as MOM4 needs thickness (and density)-weighted tendencies (multiplication by thickness)
           ! that can be skipped here.
           do k=1,grid_kmt(i,j)
              t_prog(biotic(n)%inds(ivar))%th_tendency(i,j,k) = t_prog(biotic(n)%inds(ivar))%th_tendency(i,j,k) + &
                  grid_tmask(i,j,k)*(biotic(n)%adv(k+1,ivar)-biotic(n)%adv(k,ivar)) !*Dens%rho(i,j,k,taum1)
           end do  !} k
        end if

      end do  !} ivar
    end do  !} i
  end do  !} j
enddo  !} n

call mpp_clock_end(id_clock_fabm_vertmov)

return

end subroutine  ocean_fabm_source  !}
! </SUBROUTINE> NAME="ocean_fabm_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM4 basecode
! </DESCRIPTION>
!

subroutine ocean_fabm_tracer  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

use mpp_mod, only : mpp_sum

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_tracer'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i
integer :: j
integer :: n

return

end subroutine  ocean_fabm_tracer  !}
! </SUBROUTINE> NAME="ocean_fabm_tracer"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_flux_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the ocean-atmosphere gas fluxes
! </DESCRIPTION>

subroutine ocean_fabm_flux_init  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_flux_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer                                                 :: n
character(len=fm_field_name_len)                        :: name
character(len=fm_path_name_len)                         :: path_to_names
character(len=fm_field_name_len+1)                      :: suffix
character(len=256)                                      :: caller_str

  integer :: stdoutunit
  stdoutunit=stdout()

!
!       First, perform some initialization if this module has not been
!       initialized because the normal initialization routine will
!       not have been called as part of the normal ocean model
!       initialization if this is an Atmosphere pe of a coupled
!       model running in concurrent mode
!

if (.not. module_initialized) then  !{

!
!       Initialize the package
!

  package_index = otpm_set_tracer_package(package_name,            &
       restart_file = default_restart_file,                        &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')

!
!       Check whether to use this package
!

  path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
  instances = fm_get_length(path_to_names)
  if (instances < 0) then  !{
    call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
  endif  !}

!
!       Check some things
!

  write (stdoutunit,*)
  if (instances == 0) then  !{
    write (stdoutunit,*) trim(note_header), ' No instances'
    do_ocean_fabm = .false.
  else  !}{
    if (instances == 1) then  !{
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instance'
    else  !}{
      write (stdoutunit,*) trim(note_header), ' ', instances, ' instances'
    endif  !}
    do_ocean_fabm = .true.
  endif  !}

  module_initialized = .true.

endif  !}

!
!       Return if we don't want to use this package
!

if (.not. do_ocean_fabm) then  !{
  return
endif  !}

if (.not. allocated(biotic)) then  !{

!
!       allocate storage for biotic array
!

  allocate ( biotic(instances) )

!
!       loop over the names, saving them into the biotic array
!

  do n = 1, instances  !{

    if (fm_get_value(path_to_names, name, index = n)) then  !{
      biotic(n)%name = name
    else  !}{
      write (name,*) n
      call mpp_error(FATAL, trim(error_header) //        &
           'Bad field name for index ' // trim(name))
    endif  !}

  enddo  !}

endif  !}

!
!       Set up the ocean-atmosphere gas flux fields
!

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

do n = 1, instances  !{
!  name = biotic(n)%name
!  if (name(1:1) == '_') then  !{
!    suffix = ' '
!  else  !}{
!    suffix = '_' // name
!  endif  !}
!
!!
!!       Coupler fluxes
!!
!
!  biotic(n)%ind_co2_flux = aof_set_coupler_flux('co2_flux' // suffix,                           &
!       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
!       mol_wt = WTMCO2, param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
!       ice_restart_file = default_ice_restart_file,                                             &
!       ocean_restart_file = default_ocean_restart_file,                                         &
!       caller = caller_str)
!
!  biotic(n)%ind_o2_flux = aof_set_coupler_flux('o2_flux' // suffix,                             &
!       flux_type = 'air_sea_gas_flux', implementation = 'ocmip2',                               &
!       mol_wt = WTMO2, param = (/ 9.36e-07, 9.7561e-06 /),                                                      &
!       ice_restart_file = default_ice_restart_file,                                             &
!       ocean_restart_file = default_ocean_restart_file,                                         &
!       caller = caller_str)
!
!!
!!       Coupler fields
!!

enddo  !} n

return

end subroutine  ocean_fabm_flux_init  !}
!</SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_init_sfc">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
!
!       Note: this subroutine should be merged into ocean_fabm_start
! </DESCRIPTION>

subroutine ocean_fabm_init_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,   &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     Ocean_fields, T_prog, rho, taum1, model_time, grid_tmask)  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: iec_bnd
integer, intent(in)                                     :: jsc_bnd
integer, intent(in)                                     :: jec_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), dimension(:), intent(in)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_init_sfc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer         :: i
integer         :: j
integer         :: n
integer         :: ind
real            :: epsln=1.0e-30

do n = 1, instances  !{
  allocate(biotic(n)%sf_fluxes(isc:iec,jsc:jec,1:size(biotic(n)%model%info%state_variables)))
  biotic(n)%sf_fluxes = 0.0
!
!!
!!       CO2 flux
!!
!
!  ind = biotic(n)%ind_co2_flux
!  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),    &
!                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then  !{
!
!    call ocmip2_co2calc(isd, ied, jsd, jed, isc, iec, jsc, jec,         &
!         grid_tmask(isd:ied,jsd:jed,1),                                 &
!         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                &
!         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                 &
!         t_prog(biotic(n)%ind_dic)%field(isd:ied,jsd:jed,1,taum1),      &
!         t_prog(biotic(n)%ind_alk)%field(isd:ied,jsd:jed,1,taum1),      &
!         t_prog(biotic(n)%ind_po4)%field(isd:ied,jsd:jed,1,taum1),      &
!         biotic(n)%sio2,                                                &
!         htotal_scale_lo, htotal_scale_hi, biotic(n)%htotal,            &
!         co2star = biotic(n)%csurf, alpha = biotic(n)%alpha,            &
!         pco2surf = biotic(n)%pco2surf)
!
!!
!!---------------------------------------------------------------------
!!  Compute the Schmidt number of CO2 in seawater using the
!!  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
!!  7373-7382).
!!---------------------------------------------------------------------
!!
!
!    do j = jsc, jec  !{
!      do i = isc, iec  !{
!        biotic(n)%sc_co2(i,j) =                                                 &
!             biotic(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *          &
!             (biotic(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *         &
!              (biotic(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *        &
!               biotic(n)%sc_co2_3)) * grid_tmask(i,j,1)
!        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_co2(i,j) + epsln))
!      enddo  !} i
!    enddo  !} j
!    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =   &
!         biotic(n)%alpha(isc:iec,jsc:jec) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)
!    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =   &
!         biotic(n)%csurf(isc:iec,jsc:jec) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)
!
!  endif  !}
!
!!
!!       O2 flux
!!
!
!  ind = biotic(n)%ind_o2_flux
!  if (.not. field_exist('INPUT/'//trim(Ocean_fields%bc(ind)%ocean_restart_file),    &
!                        Ocean_fields%bc(ind)%field(ind_alpha)%name)) then  !{
!
!!
!!---------------------------------------------------------------------
!!  Compute the oxygen saturation concentration at 1 atm total
!!  pressure in mol/m^3 given the temperature (t, in deg C) and
!!  the salinity (s, in permil)
!!
!!  From Garcia and Gordon (1992), Limnology and Oceonography.
!!  The formula used is from page 1310, eq (8).
!!
!!  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
!!  *** It shouldn't be there.                                ***
!!
!!  o2_saturation is defined between T(freezing) <= T <= 40 deg C and
!!                                   0 permil <= S <= 42 permil
!!
!! check value: T = 10 deg C, S = 35 permil,
!!              o2_saturation = 0.282015 mol/m^3
!!---------------------------------------------------------------------
!!
!
!    do j = jsc, jec  !{
!      do i = isc, iec  !{
!        tt(i) = 298.15 - t_prog(indtemp)%field(i,j,1,taum1)
!        tk(i) = 273.15 + t_prog(indtemp)%field(i,j,1,taum1)
!        ts(i) = log(tt(i) / tk(i))
!        ts2(i) = ts(i) * ts(i)
!        ts3(i) = ts2(i) * ts(i)
!        ts4(i) = ts3(i) * ts(i)
!        ts5(i) = ts4(i) * ts(i)
!        o2_saturation(i,j) =                                            &
!             exp(a_0 + a_1*ts(i) + a_2*ts2(i) +                         &
!                 a_3*ts3(i) + a_4*ts4(i) + a_5*ts5(i) +                 &
!                 t_prog(indsal)%field(i,j,1,taum1) *                    &
!                 (b_0 + b_1*ts(i) + b_2*ts2(i) + b_3*ts3(i) +           &
!                  c_0*t_prog(indsal)%field(i,j,1,taum1)))
!      enddo  !} i
!    enddo  !} j
!
!!
!!       convert from ml/l to mol/m^3
!!
!
!    do j = jsc, jec  !{
!      do i = isc, iec  !{
!        o2_saturation(i,j) = o2_saturation(i,j) * (1000.0/22391.6)
!      enddo  !} i
!    enddo  !} j
!
!!
!!---------------------------------------------------------------------
!!  Compute the Schmidt number of O2 in seawater using the
!!  formulation proposed by Keeling et al. (1998, Global Biogeochem.
!!  Cycles, 12, 141-163).
!!---------------------------------------------------------------------
!!
!
!    do j = jsc, jec  !{
!      do i = isc, iec  !{
!        biotic(n)%sc_o2(i,j) =                                                  &
!             biotic(n)%sc_o2_0 + t_prog(indtemp)%field(i,j,1,taum1) *           &
!             (biotic(n)%sc_o2_1 + t_prog(indtemp)%field(i,j,1,taum1) *          &
!              (biotic(n)%sc_o2_2 + t_prog(indtemp)%field(i,j,1,taum1) *         &
!               biotic(n)%sc_o2_3)) * grid_tmask(i,j,1)
!        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_o2(i,j) + epsln))
!      enddo  !} i
!    enddo  !} j
!    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         o2_saturation(isc:iec,jsc:jec) * sc_no_term(isc:iec,jsc:jec)
!    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         t_prog(biotic(n)%ind_o2)%field(isc:iec,jsc:jec,1,taum1) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)
!
!  endif  !}
enddo  !} n

return

end subroutine ocean_fabm_init_sfc  !}
! </SUBROUTINE> NAME="ocean_fabm_init_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_zero_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_fabm_zero_sfc(Ocean_fields)  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

type(coupler_2d_bc_type), intent(inout) :: Ocean_fields

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_zero_sfc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer         :: n

do n = 1, instances  !{

   biotic(n)%sf_fluxes = 0.0

!  ind = biotic(n)%ind_co2_flux
!
!  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
!  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0
!
!  ind = biotic(n)%ind_o2_flux
!
!  Ocean_fields%bc(ind)%field(ind_alpha)%values = 0.0
!  Ocean_fields%bc(ind)%field(ind_csurf)%values = 0.0

enddo  !} n

return

end subroutine ocean_fabm_zero_sfc  !}
! </SUBROUTINE> NAME="ocean_fabm_zero_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_sum_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_fabm_sum_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     Ocean_fields, T_prog, rho, taum1, model_time, grid_tmask)  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: iec_bnd
integer, intent(in)                                     :: jsc_bnd
integer, intent(in)                                     :: jec_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
type(ocean_prog_tracer_type), intent(in), dimension(:)  :: T_prog
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho
integer, intent(in)                                     :: taum1
type(time_type), intent(in)                             :: model_time
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_sum_sfc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer         :: i
integer         :: j
integer         :: n
integer         :: ivar
integer         :: ind
real            :: epsln=1.0e-30

call mpp_clock_begin(id_clock_fabm_sbc)

call get_external_fields(isc, iec, jsc, jec, nk, model_time)

do n = 1, instances  !{

  ! Set pointers to environmental variables.
  if (fabm_is_variable_used(biotic(n)%id_temp)) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_temp,t_prog(indtemp)%field(isc:iec,jsc:jec,:,taum1))
  if (fabm_is_variable_used(biotic(n)%id_salt)) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_salt,t_prog(indsal )%field(isc:iec,jsc:jec,:,taum1))
  if (fabm_is_variable_used(biotic(n)%id_dens)) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_dens,rho                  (isc:iec,jsc:jec,:,taum1))

  ! Set pointers to biotic variables.
  do ivar=1,size(biotic(n)%model%info%state_variables)
    call fabm_link_bulk_state_data(biotic(n)%model,ivar,t_prog(biotic(n)%inds(ivar))%field(isc:iec,jsc:jec,:,taum1))
  end do

  call fabm_check_ready(biotic(n)%model)

   ! Set surface fluxes for biota
   do j = jsc, jec  !{
    biotic(n)%work_dy = 0.d0
    call fabm_get_surface_exchange(biotic(n)%model,1,iec-isc+1,j-jsc+1,1,biotic(n)%work_dy(:,:))
    do ivar=1,size(biotic(n)%model%info%state_variables)
      biotic(n)%sf_fluxes(isc:iec,j,ivar) = biotic(n)%sf_fluxes(isc:iec,j,ivar) + biotic(n)%work_dy(isc:iec,ivar)*rho(isc:iec,j,1,taum1)
    end do
   enddo  !} j

!    ind = biotic(n)%ind_co2_flux
!
!    call ocmip2_co2calc(isd, ied, jsd, jed, isc, iec, jsc, jec,         &
!         grid_tmask(isd:ied,jsd:jed,1),                                 &
!         t_prog(indtemp)%field(isd:ied,jsd:jed,1,taum1),                &
!         t_prog(indsal)%field(isd:ied,jsd:jed,1,taum1),                 &
!         t_prog(biotic(n)%ind_dic)%field(isd:ied,jsd:jed,1,taum1),      &
!         t_prog(biotic(n)%ind_alk)%field(isd:ied,jsd:jed,1,taum1),      &
!         t_prog(biotic(n)%ind_po4)%field(isd:ied,jsd:jed,1,taum1),      &
!         biotic(n)%sio2,                                                &
!         htotal_scale_lo, htotal_scale_hi, biotic(n)%htotal,            &
!         co2star = biotic(n)%csurf, alpha = biotic(n)%alpha,            &
!         pco2surf = biotic(n)%pco2surf)
!
!!
!!---------------------------------------------------------------------
!!  Compute the Schmidt number of CO2 in seawater using the
!!  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
!!  7373-7382).
!!---------------------------------------------------------------------
!!
!
!    do j = jsc, jec  !{
!      do i = isc, iec  !{
!        biotic(n)%sc_co2(i,j) =                                                 &
!             biotic(n)%sc_co2_0 + t_prog(indtemp)%field(i,j,1,taum1) *          &
!             (biotic(n)%sc_co2_1 + t_prog(indtemp)%field(i,j,1,taum1) *         &
!              (biotic(n)%sc_co2_2 + t_prog(indtemp)%field(i,j,1,taum1) *        &
!               biotic(n)%sc_co2_3)) * grid_tmask(i,j,1)
!        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_co2(i,j) + epsln))
!      enddo  !} i
!    enddo  !} j
!    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) +      &
!         biotic(n)%alpha(isc:iec,jsc:jec) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)
!    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) +      &
!         biotic(n)%csurf(isc:iec,jsc:jec) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)
!
!    ind = biotic(n)%ind_o2_flux
!
!!
!!---------------------------------------------------------------------
!!  Compute the oxygen saturation concentration at 1 atm total
!!  pressure in mol/m^3 given the temperature (t, in deg C) and
!!  the salinity (s, in permil)
!!
!!  From Garcia and Gordon (1992), Limnology and Oceonography.
!!  The formula used is from page 1310, eq (8).
!!
!!  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
!!  *** It shouldn't be there.                                ***
!!
!!  o2_saturation is defined between T(freezing) <= T <= 40 deg C and
!!                                   0 permil <= S <= 42 permil
!!
!! check value: T = 10 deg C, S = 35 permil,
!!              o2_saturation = 0.282015 mol/m^3
!!---------------------------------------------------------------------
!!
!
!    do j = jsc, jec  !{
!      do i = isc, iec  !{
!        tt(i) = 298.15 - t_prog(indtemp)%field(i,j,1,taum1)
!        tk(i) = 273.15 + t_prog(indtemp)%field(i,j,1,taum1)
!        ts(i) = log(tt(i) / tk(i))
!        ts2(i) = ts(i) * ts(i)
!        ts3(i) = ts2(i) * ts(i)
!        ts4(i) = ts3(i) * ts(i)
!        ts5(i) = ts4(i) * ts(i)
!        o2_saturation(i,j) =                                        &
!             exp(a_0 + a_1*ts(i) + a_2*ts2(i) +                     &
!                 a_3*ts3(i) + a_4*ts4(i) + a_5*ts5(i) +             &
!                 t_prog(indsal)%field(i,j,1,taum1) *                &
!                 (b_0 + b_1*ts(i) + b_2*ts2(i) + b_3*ts3(i) +       &
!                  c_0*t_prog(indsal)%field(i,j,1,taum1)))
!      enddo  !} i
!    enddo  !} j
!
!!
!!       convert from ml/l to mol/m^3
!!
!
!    do j = jsc, jec  !{
!      do i = isc, iec  !{
!        o2_saturation(i,j) = o2_saturation(i,j) * (1000.0/22391.6)
!      enddo  !} i
!    enddo  !} j
!
!!
!!---------------------------------------------------------------------
!!  Compute the Schmidt number of O2 in seawater using the
!!  formulation proposed by Keeling et al. (1998, Global Biogeochem.
!!  Cycles, 12, 141-163).
!!---------------------------------------------------------------------
!!
!
!    do j = jsc, jec  !{
!      do i = isc, iec  !{
!        biotic(n)%sc_o2(i,j) =                                                    &
!             biotic(n)%sc_o2_0 + t_prog(indtemp)%field(i,j,1,taum1) *           &
!             (biotic(n)%sc_o2_1 + t_prog(indtemp)%field(i,j,1,taum1) *          &
!              (biotic(n)%sc_o2_2 + t_prog(indtemp)%field(i,j,1,taum1) *         &
!               biotic(n)%sc_o2_3)) * grid_tmask(i,j,1)
!        sc_no_term(i,j) = sqrt(660.0 / (biotic(n)%sc_o2(i,j) + epsln))
!      enddo  !} i
!    enddo  !} j
!    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) +      &
!         o2_saturation(isc:iec,jsc:jec) * sc_no_term(isc:iec,jsc:jec)
!    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) +      &
!         t_prog(biotic(n)%ind_o2)%field(isc:iec,jsc:jec,1,taum1) * rho(isc:iec,jsc:jec,1,taum1) * sc_no_term(isc:iec,jsc:jec)

enddo  !} n

call mpp_clock_end(id_clock_fabm_sbc)

return

end subroutine ocean_fabm_sum_sfc  !}
! </SUBROUTINE> NAME="ocean_fabm_sum_sfc"


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_avg_sfc">
!
! <DESCRIPTION>
!       Sum surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_fabm_avg_sfc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,    &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd, Ocean_fields, Ocean_avg_kount, grid_tmask)  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
integer, intent(in)                                     :: isc_bnd
integer, intent(in)                                     :: iec_bnd
integer, intent(in)                                     :: jsc_bnd
integer, intent(in)                                     :: jec_bnd
type(coupler_2d_bc_type), intent(inout)                 :: Ocean_fields
integer                                                 :: Ocean_avg_kount
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_avg_sfc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

integer         :: n
integer         :: ind
real            :: divid

divid = 1./float(Ocean_avg_kount)

do n = 1, instances  !{

   biotic(n)%sf_fluxes = biotic(n)%sf_fluxes * divid
!  ind = biotic(n)%ind_co2_flux
!
!  where (grid_tmask(isc:iec,jsc:jec,1) == 1.0)  !{
!    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) * divid
!    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) * divid
!  endwhere  !}
!
!  ind = biotic(n)%ind_o2_flux
!
!  where (grid_tmask(isc:iec,jsc:jec,1) == 1.0)  !{
!    Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         Ocean_fields%bc(ind)%field(ind_alpha)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) * divid
!    Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) =           &
!         Ocean_fields%bc(ind)%field(ind_csurf)%values(isc_bnd:iec_bnd,jsc_bnd:jec_bnd) * divid
!  endwhere  !}

enddo  !} n

return

end subroutine ocean_fabm_avg_sfc  !}
! </SUBROUTINE> NAME="ocean_fabm_avg_sfc"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocean_fabm_sbc(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
     isc_bnd, iec_bnd, jsc_bnd, jec_bnd,                                        &
     T_prog, taum1, model_time, grid_tmask, ice_ocean_boundary_fluxes)  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

use coupler_types_mod, only       : coupler_2d_bc_type !, ind_flux
use mpp_mod, only                 : mpp_sum

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: nk
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
integer, intent(in)                                             :: isc_bnd
integer, intent(in)                                             :: iec_bnd
integer, intent(in)                                             :: jsc_bnd
integer, intent(in)                                             :: jec_bnd
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
integer, intent(in)                                             :: taum1
type(time_type), intent(in)                                     :: model_time
real, dimension(isd:,jsd:,:), intent(in)                        :: grid_tmask
type(coupler_2d_bc_type), intent(in)                            :: ice_ocean_boundary_fluxes

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_sbc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i
integer :: j
integer :: k
integer :: n
integer :: ivar

!
! =====================================================================
!     begin executable code
! =====================================================================
!
do n = 1, instances  !{

!  ! Set pointers to environmental variables.
!  if (biotic(n)%id_temp/=-1) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_temp,t_prog(indtemp)%field(isc:iec,jsc:jec,:,taum1))
!  if (biotic(n)%id_salt/=-1) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_salt,t_prog(indsal )%field(isc:iec,jsc:jec,:,taum1))
!  !if (biotic(n)%id_dens/=-1) call fabm_link_bulk_data(biotic(n)%model,biotic(n)%id_dens,Dens%rho             (isc:iec,jsc:jec,:,taum1))
!
!  ! Set pointers to biotic variables.
!  do ivar=1,size(biotic(n)%model%info%state_variables)
!    call fabm_link_bulk_state_data(biotic(n)%model,ivar,t_prog(biotic(n)%inds(ivar))%field(isc:iec,jsc:jec,:,taum1))
!  end do

  do ivar=1,size(biotic(n)%model%info%state_variables)
    if (biotic(n)%model%info%state_variables(ivar)%no_precipitation_dilution) &
      T_prog(biotic(n)%inds(ivar))%tpme  (isc:iec,jsc:jec) = T_prog(biotic(n)%inds(ivar))%field(isc:iec,jsc:jec,1,taum1)
    if ((.not.zero_river_concentration).and.biotic(n)%model%info%state_variables(ivar)%no_river_dilution) &
      T_prog(biotic(n)%inds(ivar))%triver(isc:iec,jsc:jec) = T_prog(biotic(n)%inds(ivar))%field(isc:iec,jsc:jec,1,taum1)
  end do

!   ! Set surface fluxes for biota
!   do j = jsc, jec  !{
!    biotic(n)%work_dy = 0.d0
!    call fabm_get_surface_exchange(biotic(n)%model,1,iec-isc+1,j-jsc+1,1,biotic(n)%work_dy(:,:))
    do ivar=1,size(biotic(n)%model%info%state_variables)
      t_prog(biotic(n)%inds(ivar))%stf(isc:iec,jsc:jec) = biotic(n)%sf_fluxes(isc:iec,jsc:jec,ivar)
    end do
!   enddo  !} j
end do !} n

return

end subroutine  ocean_fabm_sbc  !}
! </SUBROUTINE> NAME="ocean_fabm_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_sfc_end">
!
! <DESCRIPTION>
!       Initialize surface fields for flux calculations
! </DESCRIPTION>

subroutine ocean_fabm_sfc_end  !{

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_sfc_end'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!       local variables
!-----------------------------------------------------------------------
!

return

end subroutine ocean_fabm_sfc_end  !}
! </SUBROUTINE> NAME="ocean_fabm_sfc_end"

!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocean_fabm_bbc(isc, iec, jsc, jec, isd, ied, jsd, jed, T_prog, grid_kmt)  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                             :: isc
integer, intent(in)                                             :: iec
integer, intent(in)                                             :: jsc
integer, intent(in)                                             :: jec
integer, intent(in)                                             :: isd
integer, intent(in)                                             :: ied
integer, intent(in)                                             :: jsd
integer, intent(in)                                             :: jed
type(ocean_prog_tracer_type), intent(inout), dimension(:)       :: T_prog
integer, dimension(isd:,jsd:), intent(in)                       :: grid_kmt

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_bbc'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer  :: i, j

!
! =====================================================================
!     begin executable code
! =====================================================================
!

  do j = jsc, jec  !{
    do i = isc, iec  !{
      ! Calculate bottom fluxes:
      ! t_prog(biotic%ind)%btf(i,j) = xxx
    enddo  !} i
  enddo  !} j

return

end subroutine  ocean_fabm_bbc  !}
! </SUBROUTINE> NAME="ocean_fabm_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_fabm_end">
!
! <DESCRIPTION>
!     Clean up various BIOTIC quantities for this run.
! </DESCRIPTION>
!

subroutine ocean_fabm_end(isc, iec, jsc, jec, nk, isd, ied, jsd, jed,        &
     T_prog, grid_dat, grid_tmask, mpp_domain2d, rho_dzt, taup1)  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

implicit none

!
!-----------------------------------------------------------------------
!       Arguments
!-----------------------------------------------------------------------
!

integer, intent(in)                                     :: isc
integer, intent(in)                                     :: iec
integer, intent(in)                                     :: jsc
integer, intent(in)                                     :: jec
integer, intent(in)                                     :: nk
integer, intent(in)                                     :: isd
integer, intent(in)                                     :: ied
integer, intent(in)                                     :: jsd
integer, intent(in)                                     :: jed
type(ocean_prog_tracer_type), intent(in), dimension(:)  :: T_prog
integer, intent(in)                                     :: taup1
real, dimension(isd:,jsd:), intent(in)                  :: grid_dat
real, dimension(isd:,jsd:,:), intent(in)                :: grid_tmask
type(domain2d), intent(in)                              :: mpp_domain2d
real, dimension(isd:,jsd:,:,:), intent(in)              :: rho_dzt

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_fabm_end'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer :: i
integer :: j
integer :: k
integer :: n

  integer :: stdoutunit
  stdoutunit=stdout()

!
!-----------------------------------------------------------------------
!     statement functions
!-----------------------------------------------------------------------
!
!
! =====================================================================
!     begin executable code
! =====================================================================
!

!
!       integrate the total concentrations of some tracers
!       for the end of the run
!

!total_phosphate = 0.0

!
!       Use taup1 time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
!

write (stdoutunit,*) trim(note_header),                           &
     'Global integrals at end of run'

do n = 1, instances  !{
  do k = 1,nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
!        total_phosphate = total_phosphate +                     &
!             t_prog(biotic(n)%ind_po4)%field(i,j,k,taup1) *     &
!             grid_dat(i,j) * grid_tmask(i,j,k) * rho_dzt(i,j,k,taup1)
      enddo  !} i
    enddo  !} j
  enddo  !} k

!  call mpp_sum(total_phosphate)

  write (stdoutunit,*) '  Instance ', trim(biotic(n)%name)
!  write (stdoutunit,                                              &
!       '(/'' Total phosphate  = '',es19.12,'' Gmol'')')         &
!       total_phosphate * 1.0e-09
enddo  !} n


return
end subroutine  ocean_fabm_end  !}
! </SUBROUTINE> NAME="ocean_fabm_end"

   subroutine mom_driver_fatal_error(self,location,message)
     class (type_mom_driver), intent(inout) :: self
     character(len=*),        intent(in)    :: location,message

     call mpp_error(FATAL, trim(location)//': '//trim(message))
   end subroutine

   subroutine mom_driver_log_message(self,message)
     class (type_mom_driver), intent(inout) :: self
     character(len=*),        intent(in)    :: message

     write (stdout(),*) trim(message)
   end subroutine

end module  ocean_fabm_mod  !}

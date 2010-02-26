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
! MOM4 driver for the Repository of Marine Biogeochemical Models (RMBM)
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
! $Id$
!

!
!------------------------------------------------------------------
!
!       Module ocean_tracer_rmbm_mod
!
!       TODO
!
!------------------------------------------------------------------
!

module  ocean_tracer_rmbm_mod  !{

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

use diag_manager_mod,         only: send_data
use field_manager_mod,        only: fm_field_name_len, fm_path_name_len, fm_string_len
use field_manager_mod,        only: fm_get_length, fm_get_value, fm_new_value
use fms_mod,                  only: write_data
use mpp_mod,                  only: stdout, stdlog, mpp_error, mpp_sum, FATAL
use time_manager_mod,         only: get_date
use time_interp_external_mod, only: time_interp_external

use ocean_tpm_util_mod, only: otpm_set_tracer_package, otpm_set_prog_tracer, otpm_set_diag_tracer
use ocean_tpm_util_mod, only: otpm_check_for_bad_fields, otpm_set_value
use ocean_tpm_util_mod, only: otpm_get_string, otpm_get_integer_array, otpm_get_logical, otpm_get_integer, otpm_get_real
use ocean_tpm_util_mod, only: otpm_get_logical_array, otpm_get_real_array, otpm_get_string_array
use ocean_tpm_util_mod, only: otpm_start_namelist, otpm_end_namelist
use ocean_tpm_util_mod, only: domain, grid, time, dtts
use ocean_tpm_util_mod, only: isc, iec, jsc, jec, nk, isd, ied, jsd, jed 
use ocean_tpm_util_mod, only: taum1, tau, taup1 
use ocean_tpm_util_mod, only: t_prog, t_diag
use ocean_tpm_util_mod, only: indsal, indtemp
use ocean_tpm_util_mod, only: end_of_year, end_of_month

use ocean_sbc_mod,      only: use_waterflux
use ocean_types_mod,    only: ocean_thickness_type

use rmbm

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

public  :: ocean_tracer_rmbm_bbc
public  :: ocean_tracer_rmbm_end
public  :: ocean_tracer_rmbm_init
public  :: ocean_tracer_rmbm_sbc
public  :: ocean_tracer_rmbm_source
public  :: ocean_tracer_rmbm_start
public  :: ocean_tracer_rmbm_tracer

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

character(len=fm_field_name_len), parameter     :: package_name = 'ocean_tracer_rmbm'
character(len=48), parameter                    :: mod_name = 'ocean_tracer_rmbm_mod'
character(len=fm_string_len), parameter         :: default_file_in = 'INPUT/ocean_tracer_rmbm.res.nc'
character(len=fm_string_len), parameter         :: default_file_out = 'RESTART/ocean_tracer_rmbm.res.nc'

!
!----------------------------------------------------------------------
!
!       Private types
!
!----------------------------------------------------------------------
!
 

type biotic_type  !{

  type (type_model),pointer :: model
  
  character(len=fm_field_name_len)            :: name
  logical                                     :: do_virtual_flux = .false.
  integer,_ALLOCATABLE,dimension(:)           :: inds,inds_diag
  double precision,pointer,dimension(:,:,:,:) :: work_state _NULL,work_diag _NULL
  double precision,_ALLOCATABLE,dimension(:,:):: w,adv
  double precision,_ALLOCATABLE,dimension(:)  :: diag _NULL,dy _NULL
  double precision,_ALLOCATABLE,dimension(:,:):: diag1d _NULL,dy1d _NULL

end type biotic_type  !}

!
!----------------------------------------------------------------------
!
!       Public variables
!
!----------------------------------------------------------------------
!

logical, public :: do_ocean_tracer_rmbm

!
!----------------------------------------------------------------------
!
!       Private variables
!
!----------------------------------------------------------------------
!

integer                                 :: package_index

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
integer                                         :: instances
type(biotic_type), allocatable, dimension(:)    :: biotic
integer                                         :: index_irr,index_chl
double precision                                :: dt_bio

!
!-----------------------------------------------------------------------
!
!       Subroutine and function definitions
!
!-----------------------------------------------------------------------
!

contains


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_rmbm_bbc">
!
! <DESCRIPTION>
!     calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocean_tracer_rmbm_bbc  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tracer_rmbm_bbc'
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

end subroutine  ocean_tracer_rmbm_bbc  !}
! </SUBROUTINE> NAME="ocean_tracer_rmbm_bbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_rmbm_end">
!
! <DESCRIPTION>
!     Clean up various BIOTIC quantities for this run.
! </DESCRIPTION>
!

subroutine ocean_tracer_rmbm_end(Thickness)  !{ 

type(ocean_thickness_type), intent(in) :: Thickness

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tracer_rmbm_end'
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
real    :: total_salinity
real    :: total_temp

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

total_temp = 0.0
total_salinity = 0.0

!
!       Use tau time index for the start of a run, and taup1 time
!       index for the end of a run so that we are integrating the
!       same time level and should therefore get identical results
!

do k = 1, nk  !{
  do j = jsc, jec  !{
    do i = isc, iec  !{
      total_temp = total_temp +                                 &
           t_prog(indtemp)%field(i,j,k,taup1) *                 &
           grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
      total_salinity = total_salinity +                         &
           t_prog(indsal)%field(i,j,k,taup1) *                  &
           grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
    enddo  !} i
  enddo  !} j
enddo  !} k

call mpp_sum(total_temp)
call mpp_sum(total_salinity)

write (stdout(),*) trim(note_header),                           &
     'Global integrals at end of run'
write (stdout(),'(/'' Total temperature  = '',es19.12,          &
                  &'' Gdeg-C m^3'')')                           &
            total_temp * 1.0e-15
write (stdout(),'(/'' Total salinity  = '',es19.12,             &
                  &'' Gsal m^3'')')                             &
            total_salinity * 1.0e-15

do k = 1,nk  !{
 do j = jsc, jec  !{
   do i = isc, iec  !{
     ! Sum variable
     !total_phosphate = total_phosphate +                     &
     !     t_prog(biotic(n)%ind_po4)%field(i,j,k,taup1) *     &
     !     grid%dat(i,j) * grid%tmask(i,j,k) * Thickness%dht(i,j,k,tau)
   enddo  !} i
 enddo  !} j
enddo  !} k

  !call mpp_sum(total_phosphate)

  !write (stdout(),*) '  Instance ', trim(biotic(n)%name)
  !write (stdout(),                                              &
  !     '(/'' Total phosphate  = '',es19.12,'' Gmol-P m^3'')')   &
  !     total_phosphate * 1.0e-09

return
end subroutine  ocean_tracer_rmbm_end  !}
! </SUBROUTINE> NAME="ocean_tracer_rmbm_end"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_rmbm_sbc">
!
! <DESCRIPTION>
!     Calculate the surface boundary conditions
! </DESCRIPTION>
!

subroutine ocean_tracer_rmbm_sbc(robert)  !{

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!

use mpp_mod, only : mpp_sum
use time_interp_external_mod, only: time_interp_external
use time_manager_mod, only: days_in_year, days_in_month,        &
     get_date, set_date

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

real, intent(in)        :: robert       ! robert time-filter coefficient

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tracer_rmbm_sbc'
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

!
! =====================================================================
!     begin executable code
! =====================================================================
!

do n = 1, instances  !{
   if (use_waterflux) then  !{

     ! Copy current surface state to rain and river input to prevent dilution due to freshwater.
     !if (biotic(n)%do_dic_virtual_flux) then  !{
     !  do j = jsc, jec  !{
     !    do i = isc, iec  !{
     !      T_prog(biotic(n)%ind_dic)%tpme(i,j) =                 &
     !           T_prog(biotic(n)%ind_dic)%field(i,j,1,tau)
     !      T_prog(biotic(n)%ind_dic)%triver(i,j) =               &
     !           T_prog(biotic(n)%ind_dic)%field(i,j,1,tau)
     !    enddo  !} i
     !  enddo  !} j
     !endif  !}

   endif  !}

   ! Set surface fluxes for biota
   do j = jsc, jec  !{
    do i = isc, iec  !{
      !t_prog(biotic(n)%ind_dic)%stf(i,j) = kw_co2(i,j) *        &
      !      biotic(n)%csat_csurf(i,j)
      !t_prog(biotic(n)%ind_o2)%stf(i,j) = kw_o2(i,j) *          &
      !      (o2_saturation(i,j) * patm_t(i,j) -                 &
      !       t_prog(biotic(n)%ind_o2)%field(i,j,1,taum1))
    enddo  !} i
   enddo  !} j 

   !
   !---------------------------------------------------------------------
   !     add in the virtual fluxes as defined by equations (2) and (3)
   !     in the biotic HOWTO.
   !       Note: the factor of 1000 is to convert the delta salinity from
   !             model units to PSU
   !---------------------------------------------------------------------
   !

   if (.not. use_waterflux) then  !{

       if (biotic(n)%do_virtual_flux) then  !{
         do j = jsc, jec  !{
           do i = isc, iec  !{
             !biotic(n)%vstf_dic(i,j) =                             &
             !     t_prog(indsal)%stf(i,j) *                        &
             !     biotic(n)%dic_global / biotic(n)%sal_global
             !t_prog(biotic(n)%ind_dic)%stf(i,j) =                  &
             !     t_prog(biotic(n)%ind_dic)%stf(i,j) +             &
             !     biotic(n)%vstf_dic(i,j)
           enddo  !} i
         enddo  !} j
       endif  !}

   endif  !}
enddo

return

end subroutine  ocean_tracer_rmbm_sbc  !}
! </SUBROUTINE> NAME="ocean_tracer_rmbm_sbc"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_rmbm_init">
!
! <DESCRIPTION>
!       Set up any extra fields needed by the tracer packages
!
!       Save pointers to various "types", such as Grid and Domains.
! </DESCRIPTION>

subroutine ocean_tracer_rmbm_init  !{

use fms_mod, only : open_namelist_file,close_file
use diag_manager_mod, only: register_diag_field

!
!       local parameters
!

character(len=64), parameter    :: sub_name = 'ocean_tracer_rmbm_init'
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
integer                                                 :: n
integer                                                 :: models(1:0)

!integer                                                 :: n
!character(len=fm_field_name_len+1)                      :: suffix
!character(len=fm_field_name_len+3)                      :: long_suffix
!character(len=fm_field_name_len)                        :: name
character(len=fm_field_name_len+1)                      :: str
integer                                                 :: nmlunit
integer                                                 :: i
integer                                                 :: nmodels
integer, pointer, dimension(:)                          :: modelids
character(len=256)                                      :: namelist_file

type (type_model), pointer                              :: childmodel

!
!       Initialize the RMBM package
!

package_index = otpm_set_tracer_package(package_name,            &
     file_in = default_file_in, file_out = default_file_out,     &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')

!
!       Check whether to use this package
!

path_to_names = '/ocean_mod/tracer_packages/' // trim(package_name) // '/names'
instances = fm_get_length(path_to_names)
if (instances .lt. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not get number of instances')
endif  !}

!
!       Check some things
!

write (stdout(),*)
if (instances .eq. 0) then  !{
  write (stdout(),*) trim(note_header), ' No instances'
  do_ocean_tracer_rmbm = .false.
else  !}{
  if (instances .eq. 1) then  !{
    write (stdout(),*) trim(note_header), ' ', instances, ' instance'
  else  !}{
    write (stdout(),*) trim(note_header), ' ', instances, ' instances'
  endif  !}
  do_ocean_tracer_rmbm = .true.
endif  !}

!
!       Return if we don't want to use this package,
!       after changing the list back
!

if (.not. do_ocean_tracer_rmbm) then  !{
  return
endif  !}

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

if (fm_new_value('/ocean_mod/GOOD/good_namelists', package_name, append = .true.) .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) //                           &
       ' Could not add ' // trim(package_name) // ' to "good_namelists" list')
endif  !}

!
!-----------------------------------------------------------------------
!       Set up the *global* namelist (shared across all instances)
!-----------------------------------------------------------------------
!

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

!call otpm_start_namelist(package_name, '*global*', caller = caller_str, no_overwrite = .true., &
!     check = .true.)

!call otpm_end_namelist(package_name, '*global*', caller = caller_str, check = .true.)

!
!-----------------------------------------------------------------------
!       Set up the instance namelists
!-----------------------------------------------------------------------
!

do n = 1, instances  !{

!
!       create the instance namelist
!

  call otpm_start_namelist(package_name, biotic(n)%name, caller = caller_str, no_overwrite = .true., &
       check = .true.)

  call otpm_set_value('namelist_file', 'bio.nml')
  call otpm_set_value('models', models, 0)
  call otpm_set_value('dt', -1.)

  call otpm_end_namelist(package_name, biotic(n)%name, check = .true., caller = caller_str)

enddo  !} n

!
!       Check for any errors in the number of fields in the namelists for this package
!

good_list => otpm_get_string_array('/ocean_mod/GOOD/namelists/' // trim(package_name) // '/good_values',   &
     caller = trim(mod_name) // '(' // trim(sub_name) // ')')
if (associated(good_list)) then  !{
  call otpm_check_for_bad_fields('/ocean_mod/namelists/' // trim(package_name), good_list,       &
       caller = trim(mod_name) // '(' // trim(sub_name) // ')')
  deallocate(good_list)
else  !}{
  call mpp_error(FATAL,trim(error_header) // ' Empty "' // trim(package_name) // '" list')
endif  !}




do n = 1, instances  !{

  call otpm_start_namelist(package_name, biotic(n)%name, caller = caller_str)
  namelist_file = otpm_get_string ('namelist_file', caller = caller_str, scalar = .true.)
  modelids      => otpm_get_integer_array('models', caller = caller_str)
  dt_bio        = otpm_get_real('dt', caller = caller_str, scalar = .true.)

  call otpm_end_namelist(package_name, biotic(n)%name, caller = caller_str)

  name = biotic(n)%name
  if (name(1:1) .eq. '_') then  !{
    suffix = ' '
    long_suffix = ' '
  else  !}{
    suffix = '_' // name
    long_suffix = ' (' // trim(name) // ')'
  endif  !}
  
  ! Create the RMBM model tree
  biotic(n)%model => rmbm_create_model()
  do i=1,ubound(modelids,1)
    childmodel => rmbm_create_model(modelids(i),parent=biotic(n)%model)
  end do
  
  ! Allow RMBM to initialize
  nmlunit = open_namelist_file(trim(namelist_file))
  call rmbm_init(biotic(n)%model,nmlunit)
  call close_file (nmlunit)

  ! Register state variables
  allocate(biotic(n)%inds(biotic(n)%model%info%state_variable_count))
  do i=1,biotic(n)%model%info%state_variable_count
     biotic(n)%inds(i) = otpm_set_prog_tracer(                                               &
          trim(biotic(n)%model%info%variables(i)%name) // suffix,                            &
          package_name,                                                                      &
          longname = trim(biotic(n)%model%info%variables(i)%longname) // trim(long_suffix),  &
          units = trim(biotic(n)%model%info%variables(i)%units),                             &
          flux_units = trim(biotic(n)%model%info%variables(i)%units)//'/s',                   &
          caller = trim(mod_name)//'('//trim(sub_name)//')',                                 &
          min_tracer_limit = biotic(n)%model%info%variables(i)%minimum,                      &
          max_tracer_limit = biotic(n)%model%info%variables(i)%maximum)
  end do

enddo  !} n

index_chl = otpm_set_diag_tracer('chl',                                             &
 caller=trim(mod_name)//'('//trim(sub_name)//')',                                      &
 longname='chlorophyll', units='mg/m^3', name_in='chl',                           &
 conversion=1.0, offset=0.0, min_tracer=0.0, max_tracer=1.e20, &
 min_range=-10.0, max_range=100.0, const_init_tracer=.true.,const_init_value=0.0,      &
 file_in='INPUT/ocean_chl.res.nc',file_out='RESTART/ocean_chl.res.nc')

index_irr = otpm_set_diag_tracer('irr',                                             &
 caller=trim(mod_name)//'('//trim(sub_name)//')',                                      &
 longname='irradiance', units='W/m^2', name_in='irr',                           &
 conversion=1.0, offset=0.0, min_tracer=0.0, max_tracer=1.e20, &
 min_range=-10.0, max_range=100.0, const_init_tracer=.true.,const_init_value=0.0,      &
 file_in='INPUT/ocean_irr.res.nc',file_out='RESTART/ocean_irr.res.nc')

!
!       register the fields
!

do n = 1, instances  !{

  if (instances .eq. 1) then  !{
    str = ' '
  else  !}{
    str = '_' // biotic(n)%name
  endif  !}

  ! Register diagnostic variables
  allocate(biotic(n)%inds_diag(biotic(n)%model%info%diagnostic_variable_count))
  do i=1,biotic(n)%model%info%diagnostic_variable_count
     biotic(n)%inds_diag(i) = register_diag_field('ocean_model',      &
          trim(biotic(n)%model%info%diagnostic_variables(i)%name)//str, grid%tracer_axes(1:3),                       &
          Time%model_time, trim(biotic(n)%model%info%diagnostic_variables(i)%longname), &
          trim(biotic(n)%model%info%diagnostic_variables(i)%units),            &
          missing_value = -1.0e+10)
  end do

enddo  !} n

if (index_irr.eq.-1 .or. index_chl.eq.-1) then
    call mpp_error(FATAL,trim(error_header) // 'index of chlorophyll and/or irradiance equals -1.')
end if

return

99 call mpp_error(FATAL,trim(error_header) // ' Unable to open RMBM namelist file "' // trim(namelist_file) // '"')

end subroutine ocean_tracer_rmbm_init  !}
! </SUBROUTINE> NAME="ocean_tracer_rmbm_init"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_rmbm_source">
!
! <DESCRIPTION>
!     compute the source terms for the BIOTICs, including boundary
!     conditions (not done in setvbc, to minimize number
!     of hooks required in MOM base code)
! </DESCRIPTION>
!

subroutine ocean_tracer_rmbm_source(Thickness)  !{

type(ocean_thickness_type), intent(in) :: Thickness

!
!-----------------------------------------------------------------------
!     modules (have to come first)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tracer_rmbm_source'
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
integer :: it,biotic_split
double precision :: dtsb
logical :: used,valid

!
! =====================================================================
!     begin executable code
! =====================================================================
!

! Set a reasonable chlorophyll concentration
T_diag(index_chl)%field(isc:iec,jsc:jec,:) = 0.08

! Calculate biological time step
biotic_split = max(1,nint(dtts/dt_bio))
dtsb = dtts/biotic_split

do n = 1, instances  !{
  ! Copy current bio state to temporary storage (we will use this temporary array to
  ! update state at the internal time step).
  do ivar=1,biotic(n)%model%info%state_variable_count
     biotic(n)%work_state(isc:iec,jsc:jec,:,ivar) = t_prog(biotic(n)%inds(ivar))%field(isc:iec,jsc:jec,:,taum1)
  end do
  
  ! Repair bio state at the start of the time step
  do k = 1, nk  !{
    do j = jsc, jec  !{
      do i = isc, iec  !{
        valid = rmbm_check_state(biotic(n)%model,i-isc+1,j-jsc+1,k,.true.)
      end do
    end do
  end do

  ! Set array with diagnostic variables to zero, becasue values for land points will not be set.
  biotic(n)%work_diag(isc:iec,jsc:jec,:,:) = 0.d0
end do

do it=1,biotic_split
   !
   !       Loop over multiple instances
   !
   do n = 1, instances  !{

     do k = 1, nk  !{
       do j = jsc, jec  !{
         do i = isc, iec  !{
           ! Move on if we are on a land point.
           if (grid%tmask(i,j,k).eq.0.) cycle
           
           ! Initialize derivatives to zero, because the bio model will increment/decrement values rather than set them.
           biotic(n)%dy = 0.d0
           
           ! Call bio model for current grid point.
           call rmbm_do(biotic(n)%model,i-isc+1,j-jsc+1,k,biotic(n)%model%info%state_variable_count,biotic(n)%model%info%diagnostic_variable_count,biotic(n)%dy,biotic(n)%diag)
           
           ! Update current state according to supplied temporal derivatives (Forward Euler)
           biotic(n)%work_state(i,j,k,:) = biotic(n)%work_state(i,j,k,:) + dtsb*biotic(n)%dy

           ! Transfer local diagnostics to global array.
           biotic(n)%work_diag(i,j,k,:) = biotic(n)%diag
         end do  !} i
       end do  !} j
     end do  !} k

   end do  !} n
enddo   !} t

do n = 1, instances  !{

   ! Repair bio state at the end of the time step
   do k = 1, nk  !{
     do j = jsc, jec  !{
       do i = isc, iec  !{
         valid = rmbm_check_state(biotic(n)%model,i-isc+1,j-jsc+1,k,.true.)
       end do
     end do
   end do
  
   ! Use updated state variable values to calculate tendencies.
   do ivar=1,biotic(n)%model%info%state_variable_count
      t_prog(biotic(n)%inds(ivar))%th_tendency(isc:iec,jsc:jec,:) = t_prog(biotic(n)%inds(ivar))%th_tendency(isc:iec,jsc:jec,:) &
         + (biotic(n)%work_state(isc:iec,jsc:jec,:,ivar)-t_prog(biotic(n)%inds(ivar))%field(isc:iec,jsc:jec,:,taum1))/dtts*Thickness%dht(isc:iec,jsc:jec,:,tau)
   end do

  ! Save diagnostic variables
  do ivar=1,biotic(n)%model%info%diagnostic_variable_count
      if (biotic(n)%inds_diag(ivar) .gt. 0) then
        used = send_data(biotic(n)%inds_diag(ivar),                        &
             biotic(n)%work_diag(isc:iec,jsc:jec,:,ivar),                  &
             time%model_time, rmask = grid%tmask(isc:iec,jsc:jec,:))
      endif
  end do

  ! Vertical movement is applied with a first-order upwind scheme.
  do j = jsc, jec
    do i = isc, iec
    
     !  Get sinking speed over entire column for all state variables.
     do k=1,nk
       call rmbm_get_vertical_movement(biotic(n)%model,i,j,k,biotic(n)%model%info%state_variable_count,biotic(n)%w(k,:))
     end do
     
     ! Interpolate to sinking speed at interfaces
     do ivar=1,biotic(n)%model%info%state_variable_count
       biotic(n)%w(2:nk,ivar) = (biotic(n)%w(2:nk,ivar) + biotic(n)%w(1:nk-1,ivar))*0.5d0
     end do
     biotic(n)%w(1,                   :) = 0.0d0   ! Surface boundary condition
     biotic(n)%w(grid%kmt(i,j)+1:nk+1,:) = 0.0d0   ! Bottom boundary condition
     
     do ivar=1,biotic(n)%model%info%state_variable_count
     
        ! Get upstream state variable values at all interfaces.
        do k=2,nk
          if (biotic(n)%w(k,ivar)>0.) then
            ! floating
            biotic(n)%adv(k,ivar) = t_prog(biotic(n)%inds(ivar))%field(i,j,k,taum1)
          elseif (biotic(n)%w(k,ivar)<0.) then
            ! sinking
            biotic(n)%adv(k,ivar) = t_prog(biotic(n)%inds(ivar))%field(i,j,k-1,taum1)
          end if
        end do
        
        ! Calculate transport over all interfaces.
        biotic(n)%adv = biotic(n)%adv*biotic(n)%w
        
        ! Add transport to tracer tendencies (conservative formulation)
        ! Note: normally the transport terms should be divided by the layer thickness.
        ! However, as MOM4 needs thickness-weighted tendencies (multiplication by thickness)
        ! that can be skipped here.
        do k=1,nk
           t_prog(biotic(n)%inds(ivar))%th_tendency(i,j,k) = t_prog(biotic(n)%inds(ivar))%th_tendency(i,j,k) + &
               grid%tmask(i,j,k)*(biotic(n)%adv(k+1,ivar)-biotic(n)%adv(k,ivar))
        end do  !} k
        
      end do  !} ivar
    end do  !} i
  end do  !} j
enddo  !} n

return

end subroutine  ocean_tracer_rmbm_source  !}
! </SUBROUTINE> NAME="ocean_tracer_rmbm_source"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_rmbm_start">
!
! <DESCRIPTION>
! Initialize variables, read in namelists, calculate constants
! for a given run and allocate diagnostic arrays
! </DESCRIPTION>
!

subroutine ocean_tracer_rmbm_start  !{

!
!-----------------------------------------------------------------------
!       modules (have to come first)
!-----------------------------------------------------------------------
!

use time_manager_mod, only: days_in_year, days_in_month,        &
     get_date, set_date
use time_interp_external_mod, only: init_external_field
use diag_manager_mod, only: register_diag_field, diag_axis_init
use fms_mod, only : read_data

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'ocean_tracer_rmbm_start'
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
integer                                                 :: i,n=1

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

!
!-----------------------------------------------------------------------
!       save the *global* namelist values
!-----------------------------------------------------------------------
!

caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

!call otpm_start_namelist(package_name, '*global*', caller = caller_str)

! Request variable values

!call otpm_end_namelist(package_name, '*global*', caller = caller_str)
      
!
!-----------------------------------------------------------------------
!       read in the namelists for each instance
!-----------------------------------------------------------------------
!




!
!-----------------------------------------------------------------------
!     give info
!-----------------------------------------------------------------------
!

write(stdout(),*)
write(stdout(),*) trim(note_header), 'Tracer runs initialized'
write(stdout(),*)

do n=1,instances
   allocate(biotic(n)%work_state(isc:iec,jsc:jec,nk,biotic(n)%model%info%state_variable_count))
   allocate(biotic(n)%work_diag (isc:iec,jsc:jec,nk,biotic(n)%model%info%diagnostic_variable_count))
   allocate(biotic(n)%dy(biotic(n)%model%info%state_variable_count))
   allocate(biotic(n)%diag(biotic(n)%model%info%diagnostic_variable_count))
   allocate(biotic(n)%w  (nk+1,biotic(n)%model%info%diagnostic_variable_count))
   allocate(biotic(n)%adv(nk+1,biotic(n)%model%info%diagnostic_variable_count))

   ! Set pointers to environmental variables
   biotic(n)%model%environment%temp   => t_prog(indtemp  )%field(isc:iec,jsc:jec,:,taum1)
   biotic(n)%model%environment%salt   => t_prog(indsal   )%field(isc:iec,jsc:jec,:,taum1)
   biotic(n)%model%environment%par    => t_diag(index_irr)%field(isc:iec,jsc:jec,:)
   biotic(n)%model%environment%par_sf => t_diag(index_irr)%field(isc:iec,jsc:jec,1)

   ! Set pointers to biotic variables.
   do i=1,biotic(n)%model%info%state_variable_count
     biotic(n)%model%state(i)%data => biotic(n)%work_state(isc:iec,jsc:jec,:,i)
   end do
enddo

return

end subroutine  ocean_tracer_rmbm_start  !}
! </SUBROUTINE> NAME="ocean_tracer_rmbm_start"


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_rmbm_tracer">
!
! <DESCRIPTION>
!     Perform things that should be done in tracer, but are done here
! in order to minimize the number of hooks necessary in the MOM4 basecode
! </DESCRIPTION>
!

subroutine ocean_tracer_rmbm_tracer  !{

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

character(len=64), parameter    :: sub_name = 'ocean_tracer_rmbm_tracer'
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

end subroutine  ocean_tracer_rmbm_tracer  !}
! </SUBROUTINE> NAME="ocean_tracer_rmbm_tracer"

end module  ocean_tracer_rmbm_mod  !}

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
!                  MSI-ERGOM (Baltic Sea Ecosystem Model)                       
!                                                                            
!  Coded by Gennadi Lessin (PML) based upon code provided by Thomas Neumann
!  (IOW) and updates the older version of ERGOM included in GOTM-BIO.
!  This older version is still provided in FABM under the name gotm/ergom for           
!  historical/educational purposes.                                           
!                                                                            
!  A detailed description of the original model is given                     
!  in T. Neumann, W. Fennel and C. Kremp, 2002. Experimental simulations     
!  with an ecosystem model of the Baltic Sea: a nutrient load reduction      
!  experiment. Global Biogeochemical Cycles, 16 (3).                         
!  http://dx.doi.org/10.1029/2001GB001450.                                   
!                                                                            
!  The present version adds oxygen-dependent phosphorus dynamics     
!  between sediment and water and the effect of bio-resuspension, as         
!  described in T. Neumann and G. Schernewski, 2008. Eutrophication in the   
!  Baltic Sea and shifts in nitrogen fixation analyzed with a 3D ecosystem   
!  model, J. Mar. Sys., 74 (1.2), pp. 592-602. 
!
!  Revision history:
!  September 2015, by Dennis Trolle (AU):
!  Implemented a switch for choosing between fresh and marine (default) environments.
!  If "fresh" is selected, oxygen debt (negative O2) due to H2S production is disabled.
!  Added a sediment burial process, and a range of additional diagnostic variables to output, 
!  incl. chlorophyll a, oxygen and nutrients in mass concentration units. 
!  Updated yaml input file with new entries (e.g., sediment burial rate, and phytoplankton
!  carbon:chla ratios for individual phyto groups)
!  May 2016, by Dennis Trolle (AU):
!  Added option to switch on or off n-fixation by cyanobacteria
!  Added settling of diatoms to bottom sediments, where diatoms are converted to fluff once settled
!
! !INTERFACE:
   MODULE selma_phytoplankton
!
! !DESCRIPTION:
!
! !USE:
   use fabm_types

   implicit none

   private
!
! !PUBLIC_DERIVED_TYPES:
  type,extends(type_base_model),public :: type_selma_phytoplankton
      type (type_state_variable_id) :: id_c
      type (type_state_variable_id) :: id_aa,id_nn,id_po,id_o2,id_dd,id_dic
      type (type_bottom_state_variable_id) :: id_fl
      type (type_dependency_id) :: id_par
      type (type_dependency_id) :: id_temp
      type (type_bottom_dependency_id) :: id_taub
      type (type_diagnostic_variable_id) :: id_chla
      type (type_diagnostic_variable_id) :: id_GPP
      type (type_diagnostic_variable_id) :: id_NPP

      real(rk) :: imin
      real(rk) :: alpha
      logical  :: nitrogen_fixation
      real(rk) :: rfr, rfc
      real(rk) :: r0
      real(rk) :: tll
      integer :: tlim
      real(rk) :: nb
      real(rk) :: deltao
      real(rk) :: Yc
      real(rk) :: sedrate
      real(rk) :: tau_crit
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
   end type
!EOP
!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the selma/phytoplankton model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!   Here, parameter values are read and the variables exported by the model are registered with FABM
!
! !INPUT PARAMETERS:
   class(type_selma_phytoplankton),intent(inout),target :: self
   integer,               intent(in)           :: configunit

   real(rk) :: c0, wz, kc
   real(rk),parameter :: secs_per_day = 86400._rk

   call self%get_parameter(c0,         'c0',    'mmol N/m3',   'background concentration',            default=0._rk)
   call self%get_parameter(self%rfr,   'rfr',   'mol P/mol N', 'phosphorus : nitrogen ratio',         default=0.0625_rk)
   call self%get_parameter(self%rfc,   'rfc',   'mol C/mol N', 'carbon : nitrogen ratio',             default=6.625_rk)
   call self%get_parameter(self%imin,  'imin',  'W/m2',        'minimal optimal light radiation',     default=50._rk)
   call self%get_parameter(self%alpha, 'alpha', 'mmol N/m3',   'half-saturation for nutrient uptake', default=0.25_rk)
   call self%get_parameter(self%r0,    'r0',    '1/d',         'maximum growth rate',                 default=1.3_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%nitrogen_fixation,    'nitrogen_fixation', '', 'whether nitrogen fixation is used to acquire nitrogen', default=.false.)
   call self%get_parameter(self%tlim,  'tlim',  '',            'temperature limitation of growth (0: none, 1: flagellate-style, 2: cyanobacteria-style)', default=0)
   select case (self%tlim)
   case (1)
      call self%get_parameter(self%tll, 'tll', 'degrees C^2', 'half-saturation temperature, squared', default=100.0_rk)
   case (2)
      call self%get_parameter(self%tll, 'tll', 'degrees C', 'lower temperature limit', default=13.5_rk)
   end select
   call self%get_parameter(self%nb,      'nb',      '1/d', 'excretion rate', default=0.01_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%deltao,  'deltao',  '1/d', 'mortality rate', default=0.02_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%Yc,      'Yc',      'mmol C/mg Chl a', 'carbon : chlorophyll a ratio', default=6.25_rk)
   call self%get_parameter(wz,           'wz',      'm/d',  'vertical velocity (positive: upwards/floating, negative: downwards/sinking)', default=0.0_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(kc,           'kc',      'm2/mmol N', 'specific light attenuation')
   !call self%get_parameter(self%sll,'sll','PSU', 'lower salinity limit', default=1.0_rk)
   !call self%get_parameter(self%sul,'sul','PSU', 'upper salinity limit', default=10.0_rk)
   call self%get_parameter(self%sedrate, 'sedrate', 'm/d', 'sedimentation rate', default=0.0_rk, scale_factor=1.0_rk/secs_per_day)

   call self%register_state_variable(self%id_c, 'c', 'mmol N/m3', 'concentration', minimum=0.0_rk, background_value=c0, vertical_movement=wz)
   call self%register_state_dependency(self%id_aa, 'aa', 'mmol N/m3', 'ammonium')
   call self%register_state_dependency(self%id_nn, 'nn', 'mmol N/m3', 'nitrate')
   call self%register_state_dependency(self%id_o2, 'o2', 'mmol O2/m3','oxygen')
   call self%register_state_dependency(self%id_po, 'po', 'mmol P/m3', 'phosphate')
   call self%register_state_dependency(self%id_dd, 'dd', 'mmol N/m3', 'detritus')
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_par,  standard_variables%downwelling_photosynthetic_radiative_flux)
   if (self%sedrate>0.0_rk) then
      call self%get_parameter(self%tau_crit,'tau_crit','N/m2', 'critical shear stress', default=0.07_rk)
      call self%register_state_dependency(self%id_fl, 'fl', 'mmol N/m2', 'fluff')
      call self%register_dependency(self%id_taub, standard_variables%bottom_stress)
   end if
   call self%register_state_dependency(self%id_dic,standard_variables%mole_concentration_of_dissolved_inorganic_carbon, required=.false.)

   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_c)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, self%rfr)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_c, self%rfc)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_c, kc, include_background=.true.)

   call self%register_diagnostic_variable(self%id_chla, 'chla','mg chl a/m3', 'chlorophyll concentration')
   call self%register_diagnostic_variable(self%id_GPP,  'GPP', 'mmol/m3/d',   'gross primary production')
   call self%register_diagnostic_variable(self%id_NPP,  'NPP', 'mmol/m3/d',   'net primary production')

!  we create an aggregate variable for chlA
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_chlorophyll',units="mg/m^3",aggregate_variable=.true.),self%id_chla,scale_factor=1._rk)

   end subroutine initialize
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of selma model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class(type_selma_phytoplankton), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk) :: par, temp
   real(rk) :: c, cg
   real(rk) :: aa, nn, po, o2
   real(rk) :: iopt, ppi
   real(rk) :: ntemp, ptemp, tempq
   real(rk) :: nlim, plim, r
   real(rk) :: lpn, lpd
   real(rk),parameter :: epsilon = 0.00000000001_rk
   real(rk),parameter :: secs_per_day = 86400._rk

   ! Enter spatial_loops (if any)
   _LOOP_BEGIN_

      _GET_(self%id_c,c)                  ! own concentration
      _GET_WITH_BACKGROUND_(self%id_c,cg) ! own concentration including background

      _GET_(self%id_aa,aa) ! ammonium
      _GET_(self%id_nn,nn) ! nitrate
      _GET_(self%id_po,po) ! phosphate
      _GET_(self%id_o2,o2) ! oxygen

      _GET_(self%id_par,par)   ! local photosynthetically ative radiaiton (W/m2)
      _GET_(self%id_temp,temp) ! temperature (degrees Celsius)

      ! Light acclimation formulation based on surface light intensity
      iopt = max(0.5_rk*par,self%imin)
      ppi = par/iopt*exp(1.0_rk-par/iopt)

      ! Nitrogen limitation (type III fucntional response)
      if (self%nitrogen_fixation) then
         nlim = 1.0_rk
      else
         ntemp = (nn + aa)**2
         nlim = ntemp / (self%alpha * self%alpha + ntemp) ! MiMe eq. for IN
      end if

      ! Phosphorus limitation (type III fucntional response)
      ptemp = po**2
      plim = ptemp / (self%alpha * self%alpha * self%rfr * self%rfr + ptemp)

      r = min(nlim, plim, ppi) * self%r0

      ! Temperature limitation
      if (self%tlim == 1) then
         ! Flagellate-style: Type III [sigmoidal] functional response ranging between 1 and 2
         tempq = max(temp,0.0_rk)**2
         r = r * (1.0_rk + tempq / (self%tll + tempq))
      elseif (self%tlim == 2) then
         ! Cyanobacteria-style: 0 at infinitely low temp, 0.5 at tll, 1.0 at infinitely high temp
         r = r  / (1.0_rk + exp(self%tll - temp))
      end if

      lpn = self%nb     ! excretion rate
      lpd = self%deltao ! mortality rate

      if (o2 <= 0.0_rk) then
         ! Anoxic: no growth or respiration, higher mortality
         r = 0.0_rk
         lpn = 0.0_rk
         lpd = 10._rk * self%deltao
      end if

      if (self%nitrogen_fixation) then
         ! Nitrogen acquired from dinitrogen gas (not tracked)
         _ADD_SOURCE_(self%id_o2, + r * cg * 6.625_rk)
      else
         ! Nitrogen acquired from ammonium and nitrate pools (proportional to availability)
         _ADD_SOURCE_(self%id_aa, - r * cg * aa/(nn + aa + epsilon))
         _ADD_SOURCE_(self%id_nn, - r * cg * (1.0_rk - aa/(nn + aa + epsilon)))
         _ADD_SOURCE_(self%id_o2,   r * cg * (nn * 8.625_rk + aa * 6.625_rk)/(nn + aa + epsilon))
      end if

      _ADD_SOURCE_(self%id_o2, - 6.625_rk * (lpn * c))
      _ADD_SOURCE_(self%id_aa, + lpn * c)
      _ADD_SOURCE_(self%id_po, self%rfr * (- r * cg + lpn * c))
      _ADD_SOURCE_(self%id_c, r * cg - (lpn + lpd) * c)
      _ADD_SOURCE_(self%id_dd, lpd * c)

      if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic, self%rfc*(lpn * c - r * cg))

      _SET_DIAGNOSTIC_(self%id_chla, (c * 9.5 + 2.1)/self%Yc) ! relation between carbon and nitrogen from Hecky et al 1993. The stoichiometry of carbon, nitrogen, and phosphorus in particulate matter of lakes and oceans. Limnology and Oceanography, 38: 709-724.
      _SET_DIAGNOSTIC_(self%id_GPP, secs_per_day * r * cg)
      _SET_DIAGNOSTIC_(self%id_NPP, secs_per_day *(r * cg - lpn * c))

   ! Leave spatial loops (if any)
   _LOOP_END_

   END subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of benthic_predator model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_selma_phytoplankton),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: c
   real(rk)                   :: taub
   real(rk)                   :: ll
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (self%sedrate == 0.0_rk) return

   ! Enter spatial loops over the horizontal domain (if any).
   _BOTTOM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
 !   if (self%fluff) then
   _GET_(self%id_c,c)
   _GET_BOTTOM_(self%id_taub,taub)

   ! Resuspension-sedimentation rate are computed as in GOTM-BIO
   ! Phytoplankton is assumed to become detritus/fluff as soon as it settles to bottom sediments,
   ! and can therefore not be resuspended from benthic layer

   ll=self%sedrate*max(0.0_rk, self%tau_crit-taub)/self%tau_crit

   ! Sediment resuspension, detritus settling, diatom settling, bio-resuspension, mineralization and burial
   _ADD_BOTTOM_SOURCE_(self%id_fl,+ ll * c)
   _ADD_BOTTOM_FLUX_(self%id_c, -ll * c)

   ! Leave spatial loops over the horizontal domain (if any).
   _BOTTOM_LOOP_END_

   end subroutine do_bottom
!EOC

!-----------------------------------------------------------------------

  END MODULE selma_phytoplankton

!-----------------------------------------------------------------------

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!  SELMA is yet another re-write of ERGOM/MSI-ERGOM
!  SELMA has the same functionality as MSI-ERGOM but uses modular
!  components encouraged when developing FABM models.
!  SELMA was developed in the PROGNOS project http://http://prognoswater.org/ - 
!  and the re-write was done by Jorn Bruggeman.
!  Renaming ERGOM to SELMA was done in agreement with Thomas Neumann - the
!  original author of ERGOM.
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
! !MODULE: selma
!
! !INTERFACE:
   MODULE selma
!
! !DESCRIPTION:
!
! !USE:
   use fabm_types

   implicit none

   private
!
! !PUBLIC_DERIVED_TYPES:
   type,extends(type_base_model),public :: type_selma
      ! Variable identifiers
      type (type_state_variable_id)              :: id_dd,id_aa,id_nn,id_po,id_o2,id_pw,id_dic
      type (type_bottom_state_variable_id)       :: id_fl,id_pb
      type (type_dependency_id)                  :: id_temp,id_salt
      type (type_bottom_dependency_id)           :: id_taub
      type (type_surface_dependency_id)          :: id_wind
      type (type_diagnostic_variable_id)         :: id_DNP,id_NO3,id_NH4,id_PO4,id_O2_mg
      type (type_bottom_diagnostic_variable_id)  :: id_DNB,id_SBR,id_PBR
      type (type_surface_diagnostic_variable_id) :: id_OFL

      ! Model parameters
      real(rk) :: nb,deltao,nue,sigma_b,dn,dn_sed
      real(rk) :: q10_rec,ade_r0,alphaade,q10_recs
      real(rk) :: rfr,rfc,sedrate,erorate,sedratepo4,eroratepo4,po4ret
      real(rk) :: fl_burialrate,pburialrate,pliberationrate,ipo4th,maxsed,br0,fds,pvel,tau_crit
      integer  :: newflux
      character(len=16) :: env_type ! identifier for setting the environment to "marine" or "fresh" (freshwater/lake) 
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the selma model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!   Here, parameter values are read read and the variables exported by the model are registered with FABM
!
! !INPUT PARAMETERS:
   class (type_selma), intent(inout), target :: self
   integer,            intent(in)            :: configunit
   real(rk),parameter :: secs_per_day = 86400._rk
   real(rk) :: wdz,wpo4,kc
   
   call self%get_parameter(self%env_type, 'env_type', 'Define environment type, either fresh or marine', default='marine') 
   call self%get_parameter(wdz,     'wdz',   'm/d',  'vertical velocity of detritus (positive: upwards/floating, negative: downwards/sinking)', default=-4.5_rk)
   call self%get_parameter(wpo4,    'wpo4',  'm/d',  'vertical velocity of suspended P-Fe (positive: upwards/floating, negative: downwards/sinking)', default=-1.0_rk)
   call self%get_parameter(self%dn,      'dn',      '1/d', 'detritus mineralization rate', default=0.003_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%dn_sed,  'dn_sed',  '1/d', 'sediment mineralization rate', default=0.002_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(kc,      'kc',      'm2/mmol N', 'specific light attenuation of detritus')
   call self%get_parameter(self%q10_rec, 'q10_rec', '1/K', 'temperature dependence of detritus remineralization', default=0.15_rk)
   call self%get_parameter(self%ade_r0,  'ade_r0',  '1/d', 'maximum chemoautolithotrophic denitrification rate', default=0.1_rk)
   self%ade_r0 = self%ade_r0/secs_per_day
   call self%get_parameter(self%alphaade,'alphaade','mmol N/m3', 'half-saturation constant for chemoautolithotrophic denitrification', default=1.0_rk)
   call self%get_parameter(self%q10_recs,'q10_recs','1/K', 'temperature dependence of sediment remineralization', default=0.175_rk)
   call self%get_parameter(self%rfr,     'rfr',     '-', 'phosphorus : nitrogen ratio of detritus and sediment', default=0.0625_rk)
   call self%get_parameter(self%rfc,     'rfc',     '-', 'carbon : nitrogen ratio of detritus and sediment', default=6.625_rk)
   call self%get_parameter(self%tau_crit,'tau_crit','N/m2', 'critical shear stress', default=0.07_rk)
   call self%get_parameter(self%sedrate, 'sedrate', 'm/d', 'detritus sedimentation rate', default=2.25_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%erorate, 'erorate', '1/d', 'sediment erosion rate', default=6._rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%sedratepo4,'sedratepo4','m/d', 'P-Fe sedimentation rate', default=0.5_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%eroratepo4,'eroratepo4','1/d', 'P-Fe erosion rate', default=6._rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%po4ret,  'po4ret',  '-', 'phosphate retention rate, oxic sediments', default=0.18_rk)
   call self%get_parameter(self%pburialrate,'pburialrate','1/d', 'phosphate burial rate', default=0.007_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%fl_burialrate,'fl_burialrate','1/d', 'sediment burial rate', default=0.001_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%pliberationrate,'pliberationrate','1/d', 'phosphate liberation rate, anoxic sediments', default=0.1_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%ipo4th,  'ipo4th',  'mmol P/m2', 'maximum phosphorus density available for burial', default=100._rk)
   call self%get_parameter(self%maxsed,  'maxsed',  'mmol N/m2', 'maximum active sediment density', default=1000._rk)
   call self%get_parameter(self%br0,     'br0',     '1/d', 'bioresuspension rate', default=0.03_rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%fds,     'fds',     '-', 'fraction of sediment remineralization fueled by denitrification', default=0.7_rk)
   call self%get_parameter(self%pvel,    'pvel',    'm/d', 'piston velocity', default=5._rk, scale_factor=1.0_rk/secs_per_day)
   call self%get_parameter(self%newflux, 'newflux', '-', 'oxygen flux type', default=2)
   
   ! Register state variables
   call self%register_state_variable(self%id_dd,'dd','mmol N/m3', 'detritus', minimum=0.0_rk,vertical_movement=wdz/secs_per_day,no_river_dilution=.true.)
   call self%register_state_variable(self%id_aa,'aa','mmol N/m3', 'ammonium', minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_nn,'nn','mmol N/m3', 'nitrate', minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_po,'po','mmol P/m3', 'phosphate', minimum=0.0_rk,no_river_dilution=.true.)
   if (self%env_type .eq. "fresh") then
      call self%register_state_variable(self%id_o2,'o2','mmol O2/m3','oxygen', minimum=0.0_rk,no_river_dilution=.true.)
   else
      call self%register_state_variable(self%id_o2,'o2','mmol O2/m3','oxygen', no_river_dilution=.true.)
   end if
   call self%register_state_variable(self%id_fl,'fl','mmol N/m2', 'fluff', minimum=0.0_rk,maximum=self%maxsed)
   call self%register_state_variable(self%id_pb,'pb','mmol P/m2', 'PFe_s', minimum=0.0_rk)
   call self%register_state_variable(self%id_pw,'pw','mmol P/m3', 'PFe_w', minimum=0.0_rk,vertical_movement=wpo4/secs_per_day,no_river_dilution=.true.)
     
   ! Register optional link to external DIC pool.
   call self%register_state_dependency(self%id_dic,standard_variables%mole_concentration_of_dissolved_inorganic_carbon, required=.false.)

   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_aa)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_nn)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_dd)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_dd, scale_factor=self%rfr)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_dd, scale_factor=self%rfc)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_po)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_fl)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_fl, scale_factor=self%rfr)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_fl, scale_factor=self%rfc)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_dd, kc)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_NO3,    'Nit',      'mg NO3N/m**3','nitrate conc in mass unit')
   call self%register_diagnostic_variable(self%id_NH4,    'Amm',      'mg NH4N/m**3','ammonium  conc in nitrogen mass unit')
   call self%register_diagnostic_variable(self%id_PO4,    'Pho',      'mg PO4P/m**3','phosphate conc in phosphorus mass unit')
   call self%register_diagnostic_variable(self%id_O2_mg,  'DO_mg',    'mg O2/m**3',  'oxygen in O2 mass unit')
   call self%register_diagnostic_variable(self%id_DNP,     'DNP',     'mmol N/m3/d','denitrification pelagic')
   call self%register_diagnostic_variable(self%id_DNB,     'DNB',     'mmol N/m2/d','denitrification benthic')
   call self%register_diagnostic_variable(self%id_SBR,     'SBR',     'mmol N/m2',  'sediment burial')
   call self%register_diagnostic_variable(self%id_PBR,     'PBR',     'mmol P/m2/d','phosphorus burial')
   call self%register_diagnostic_variable(self%id_OFL,     'OFL',     'mmol O2/m2/d','oxygen surface flux (positive when into water)')

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
   call self%register_dependency(self%id_wind, standard_variables%wind_speed)
   call self%register_dependency(self%id_taub, standard_variables%bottom_stress)

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
   class(type_selma), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
    real(rk)           :: aa,nn,po,dd,o2
    real(rk)           :: temp,ldn,nf
    real(rk)           :: otemp,ldn_N,ldn_O,ade
    real(rk),parameter :: p_molar_mass = 30.973761_rk ! molar mass of phosphorus
    real(rk),parameter :: n_molar_mass =  14.0067_rk ! molar mass of nitrogen
    real(rk),parameter :: o2_molar_mass =  31.9988_rk ! molar mass of O2
    real(rk),parameter :: epsilon = 0.00000000001_rk
    real(rk),parameter :: secs_per_day = 86400._rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial_loops (if any)
   _LOOP_BEGIN_

      ! Retrieve current (local) state variable values
      _GET_(self%id_dd,dd) ! detritus
      _GET_(self%id_aa,aa) ! ammonium
      _GET_(self%id_nn,nn) ! nitrate
      _GET_(self%id_po,po) ! phosphate
      _GET_(self%id_o2,o2) ! oxygen

      _GET_(self%id_temp,temp)

      !!!! NITRIFICATION RATE !!!!
      otemp = max(0.0_rk, o2) !for oxygen dependent process if o2<0 then o2=0
      ! Nitrification rate depends on oxygen availability and temperature
      nf = otemp / (0.01_rk + otemp) * 0.1_rk * exp (0.11_rk * temp)/secs_per_day

      ldn = self%dn * exp (self%q10_rec*temp) ! Mineralization rate depends on temperature

      if (o2 <= 0.0_rk .and. nn > 0.0_rk) then
         ldn_N = 5.3_rk * nn * nn / (0.001_rk + nn * nn)               ! Denitrification rate depends on nitrate availability
         ldn_O = 6.625_rk * (1.0_rk - nn * nn / (0.001_rk + nn * nn))  ! Oxygen loss due to denitrification
         ade = self%ade_r0 * nn * nn / (self%alphaade +  nn * nn) * nn ! ade rate nitrate dependent
      else
         ldn_N = 0.0_rk
         ldn_O = 6.625_rk  ! Negative oxygen = H2S production
         ade = 0.0_rk
      endif

      _ADD_SOURCE_(self%id_o2, - ldn_O * ldn * dd - 2.0_rk * nf * aa + ade * 0.375_rk)
      _ADD_SOURCE_(self%id_aa, ldn * dd - nf * aa)
      _ADD_SOURCE_(self%id_nn, nf * aa - ldn_N * ldn * dd - ade)
      _ADD_SOURCE_(self%id_po, self%rfr * ldn * dd)
      _ADD_SOURCE_(self%id_dd, - ldn * dd)

      if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic, self%rfc * ldn * dd)

      _SET_DIAGNOSTIC_(self%id_NO3, nn * n_molar_mass)
      _SET_DIAGNOSTIC_(self%id_NH4, aa * n_molar_mass)
      _SET_DIAGNOSTIC_(self%id_PO4, po * p_molar_mass)
      _SET_DIAGNOSTIC_(self%id_O2_mg, o2 * o2_molar_mass)
      _SET_DIAGNOSTIC_(self%id_DNP, (ldn_N * ldn *dd + ade) * secs_per_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
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
   class (type_selma),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: pwb,pb,fl,nnb,ddb,oxb,taub,temp,biores
   real(rk)                   :: llds,llsd,bpds,bpsd,recs,ldn_O,ldn_N,plib,oxlim
   real(rk)                   :: fracdenitsed,pret,pbr
   real(rk)                   :: tau_frac
   real(rk),parameter :: secs_per_day = 86400._rk
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _BOTTOM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
 !   if (self%fluff) then
   _GET_(self%id_dd,ddb)
   _GET_(self%id_nn,nnb)
   _GET_(self%id_o2,oxb)
   _GET_(self%id_pw,pwb)
   _GET_BOTTOM_(self%id_fl,fl)
   _GET_BOTTOM_(self%id_pb,pb)

   _GET_BOTTOM_(self%id_taub,taub)
   _GET_(self%id_temp,temp)

   !increased phosphorus burial
   pbr = max(pb, pb * (pb -self%ipo4th + 1.0_rk)) 

   !bio-resuspension rate
   biores = self%br0 * max(0.0_rk, oxb) * max(0.0_rk, oxb) &
            / (max(0.0_rk, oxb) * max(0.0_rk, oxb) + 0.03_rk)  

   ! Resuspension-sedimentation rate are computed as in GOTM-BIO
  
   tau_frac = (taub-self%tau_crit)/self%tau_crit
   ! if actual taub is greater than critical tau, then do resuspension
   if (taub .gt. self%tau_crit) then
      llds=0.0_rk
      llsd=self%erorate*(taub-self%tau_crit)/self%tau_crit
      bpds=0.0_rk
      bpsd=self%eroratepo4*(taub-self%tau_crit)/self%tau_crit
   else
      llds=self%sedrate*(self%tau_crit-taub)/self%tau_crit
      llsd=0.0_rk
      bpds=self%sedratepo4*(self%tau_crit-taub)/self%tau_crit
      bpsd=0.0_rk
   end if

   recs = self%dn_sed * exp(self%q10_recs * temp)

   ! temp-dependent detritus mineralization rate
   if (oxb < 0.0_rk) then
      ! 10 times lower in anoxic
      recs = recs * 0.1_rk
   endif

   ! Mineralization rates (see description of pelagic part)
   if (oxb .le. 0.0_rk .and. nnb > 0.0_rk) then
      ldn_N = 5.3_rk * nnb * nnb / (0.001_rk + nnb * nnb)
      ldn_O = 6.625_rk * (1.0_rk - nnb * nnb / (0.001_rk + nnb * nnb))
   else
      ldn_N = 0.0_rk
      ldn_O = 6.625_rk
   endif

   if (oxb > 0.0_rk) then
      fracdenitsed = self%fds     ! denitrification in sediments
      plib = 0.0_rk               ! no phosphate release
      pret = self%po4ret          ! phosphate is stored
   else
      fracdenitsed = 0.0_rk       ! no denitrification in sediments
      plib = self%pliberationrate ! phosphorus is liberated
      pret = 0.0_rk               ! no phosphate retention
   endif

   oxlim = max (0.0_rk,oxb) * max (0.0_rk,oxb) / (0.01_rk + max(0.0_rk,oxb) * max(0.0_rk,oxb))

   ! Sediment resuspension, detritus settling, diatom settling, bio-resuspension, mineralization and burial
   _ADD_BOTTOM_SOURCE_(self%id_fl,-llsd * fl + llds * ddb - biores * fl - recs * fl - fl * self%fl_burialrate * fl/self%maxsed)
   ! P-Fe resuspension, sedimentation, bio-resuspension, liberation, retention and burial
   _ADD_BOTTOM_SOURCE_(self%id_pb,-bpsd * pb + bpds * pwb -biores * pb -plib * pb + self%rfr * recs * fl * pret * oxlim - pbr * self%pburialrate * fl/self%maxsed)

   ! Denitrification in sediments
   _ADD_BOTTOM_FLUX_(self%id_nn,-ldn_N * recs * fl)
   ! Oxygen consumption due to mineralization and denitrification
   _ADD_BOTTOM_FLUX_(self%id_o2,-ldn_O * recs * fl - 2.0_rk * fracdenitsed * recs * fl)
   ! Ammonium production due to mineralization (oxic & anoxic)
   _ADD_BOTTOM_FLUX_(self%id_aa,(1.0_rk - fracdenitsed) * recs * fl)
   ! Phosphate production due to mineralization (retention if oxic) and release in anoxic
   _ADD_BOTTOM_FLUX_(self%id_po, (1.0_rk - pret * oxlim) * self%rfr * recs * fl + plib * pb)
   ! Sediment resuspension, detritus settling, bio-resuspension
   _ADD_BOTTOM_FLUX_(self%id_dd, llsd * fl -llds * ddb + biores * fl)
   ! P-Fe resuspension, settling and bio-resuspension
   _ADD_BOTTOM_FLUX_(self%id_pw, bpsd * pb -bpds * pwb + biores * pb)

   if (_AVAILABLE_(self%id_dic)) _ADD_BOTTOM_FLUX_(self%id_dic, self%rfc*recs * fl)

   ! BENTHIC DIAGNOSTIC VARIABLES
   _SET_BOTTOM_DIAGNOSTIC_(self%id_DNB,(ldn_N * recs * fl + fracdenitsed * recs * fl) * secs_per_day)
   _SET_BOTTOM_DIAGNOSTIC_(self%id_SBR,(fl * self%fl_burialrate * fl/self%maxsed) * secs_per_day)
   _SET_BOTTOM_DIAGNOSTIC_(self%id_PBR,(pbr * self%pburialrate * fl/self%maxsed) * secs_per_day)

   ! Leave spatial loops over the horizontal domain (if any).
   _BOTTOM_LOOP_END_

   end subroutine do_bottom
!EOC

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Weiss formula for the saturation oxygen (osat)
!
! !INTERFACE:
   real(rk) function osat_weiss(t,s)
!
! !DESCRIPTION:
! Weiss formula for the saturation oxygen (osat) \cite{Weiss1970}:
!
! !INPUT PARAMETERS:
!  type (type_selma), intent(in) :: self
  real(rk), intent(in)                 :: t,s
!
! !LOCAL VARIABLES:
  real(rk)                 :: tk
  real(rk),parameter       :: aa1=-173.4292_rk
  real(rk),parameter       :: aa2=249.6339_rk
  real(rk),parameter       :: a3=143.3483_rk
  real(rk),parameter       :: a4=-21.8492_rk
  real(rk),parameter       :: b1=-0.033096_rk
  real(rk),parameter       :: b2=0.014259_rk
  real(rk),parameter       :: b3=-0.001700_rk
  real(rk),parameter       :: kelvin=273.16_rk
  real(rk),parameter       :: mol_per_liter=44.661_rk

!EOP
!-----------------------------------------------------------------------
!BOC
   tk=(t+kelvin)*0.01_rk
   osat_weiss=exp(aa1+aa2/tk+a3*log(tk)+a4*tk    &
              +s*(b1+(b2+b3*tk)*tk))*mol_per_liter

   end function osat_weiss
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the selma model
!
! !INTERFACE:
  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !INPUT PARAMETERS:
  class(type_selma),intent(in) :: self

  _DECLARE_ARGUMENTS_DO_SURFACE_

  real(rk)           :: temp,wnd,salt,o2
  real(rk),parameter :: secs_per_day = 86400._rk
!
! !LOCAL VARIABLES:
  real(rk)                 :: p_vel,sc,flo2,osat
!  integer,parameter        :: newflux=2
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
   _SURFACE_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_SURFACE_(self%id_wind,wnd)

   _GET_(self%id_o2,o2)

!  Calculation of the surface oxygen flux
! Newflux=1 piston velocity depends on winds and saturation is
! calculated according to Weiss formula
   if (self%newflux .eq. 1) then
      sc=1450._rk+(1.1_rk*temp-71._rk)*temp
      if (wnd .gt. 13._rk) then
         p_vel = 5.9_rk*(5.9_rk*wnd-49.3_rk)/sqrt(sc)
      else
         if (wnd .lt. 3.6_rk) then
            p_vel = 1.003_rk*wnd/(sc)**(0.66_rk)
         else
            p_vel = 5.9_rk*(2.85_rk*wnd-9.65_rk)/sqrt(sc)
         end if
      end if
      if (p_vel .lt. 0.05_rk) then
          p_vel = 0.05_rk
      end if
      p_vel = p_vel/secs_per_day
      flo2 =p_vel*(osat_weiss(temp,salt)-o2)
      _ADD_SURFACE_FLUX_(self%id_o2,flo2)
! Newflux=2 use polynom approximation
! to order S & T ** 3 for oxygen saturation derived from
! www.helcom.fi/Monas/CombineManual2/PartB/AnnexB-8Appendix3.pdf
   elseif (self%newflux .eq. 2) then
      osat = (10.18e0_rk + ((5.306e-3_rk - 4.8725e-5_rk * temp) *temp - 0.2785e0_rk) * temp &
          + salt * ((2.2258e-3_rk + (4.39e-7_rk * temp - 4.645e-5_rk) * temp) * temp - 6.33e-2_rk)) &
          * 44.66e0_rk
      flo2 = self%pvel * (osat - o2)
      _ADD_SURFACE_FLUX_(self%id_o2,flo2)
   else
      flo2 = self%pvel * (31.25_rk * (14.603_rk - 0.40215_rk * temp) - o2)
      _ADD_SURFACE_FLUX_(self%id_o2,flo2)
   end if

   _SET_SURFACE_DIAGNOSTIC_(self%id_OFL,flo2 * secs_per_day)

   _SURFACE_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------

  END MODULE selma

!-----------------------------------------------------------------------

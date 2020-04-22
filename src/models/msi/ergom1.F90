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
! !MODULE:fabm_msi_ergom1
!
! !INTERFACE:
   MODULE fabm_msi_ergom1
!
! !DESCRIPTION:
!
! !USE:
   use fabm_types

   implicit none

   private
!
! !PUBLIC_DERIVED_TYPES:
  type,extends(type_base_model),public :: type_msi_ergom1
      ! Variable identifiers
      type (type_state_variable_id)        :: id_pp,id_ff,id_bb,id_zz,id_dd,id_aa,id_nn,id_po,id_o2,id_pw,id_dic
      type (type_bottom_state_variable_id) :: id_fl,id_pb
      type (type_dependency_id)            :: id_par,id_temp,id_salt
      type (type_horizontal_dependency_id) :: id_I_0,id_taub,id_wind
      type (type_diagnostic_variable_id)   :: id_dPAR,id_GPP,id_NCP,id_PPR,id_NPR,id_NFX,id_DNP,id_NO3,id_NH4,id_TN,id_PO4,id_TP,id_O2_mg,id_pp_chla,id_ff_chla,id_bb_chla,id_tot_chla
      type (type_bottom_diagnostic_variable_id) :: id_DNB,id_SBR,id_PBR
      type (type_surface_diagnostic_variable_id) :: id_OFL

      ! Model parameters
      real(rk) :: nb,deltao,nue,sigma_b,dn,dn_sed,sfl_po,sfl_aa,sfl_nn
      real(rk) :: rp0,rf0,rb0,cyanotll,cyanosll,cyanosul,flagtll,Yc_cyan,Yc_diat,Yc_flag
      real(rk) :: alphap,alphaf,alphab,iv,graz,toptz,zcl1,p0,f0,b0,z0
      real(rk) :: imin_di,imin,kc,q10_rec,ade_r0,alphaade,q10_recs
      real(rk) :: rfr,rfc,sedrate,erorate,sedratepo4,eroratepo4,po4ret
      real(rk) :: fl_burialrate,pburialrate,pliberationrate,ipo4th,maxsed,br0,fds,pvel,tau_crit
      integer  :: newflux
      logical  :: calc_dic=.false.
      character(len=16) :: env_type ! identifier for setting the environment to "marine" or "fresh" (freshwater/lake) 
      character(len=64) :: dic_variable      
      logical  :: nitrogen_fixation=.true. ! identifier for whether or not to allow n-fixation for cyanos "true" to use n-fix or "false" to allow n-limitation function
   contains
      procedure :: initialize
      procedure :: do
      procedure :: get_light_extinction
      procedure :: do_surface
      procedure :: do_bottom
   end type
!EOP
!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the ergom1 model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!   Here, parameter values are read and the variables exported by the model are registered with FABM
!
! !INPUT PARAMETERS:
   class(type_msi_ergom1),intent(inout),target :: self
   integer,               intent(in)           :: configunit
   real(rk),parameter :: secs_per_day = 86400._rk
   real(rk) :: wdz,wpz,wfz,wbz,wpo4
   
   call self%get_parameter(self%env_type, 'env_type', 'Define environment type, either fresh or marine', default='marine') 
   call self%get_parameter(self%calc_dic, 'calc_dic', 'Deside whether to calculate DIC', default=.false.)
   call self%get_parameter(self%dic_variable, 'dic_variable', 'Define DIC variable')
   call self%get_parameter(self%nitrogen_fixation, 'nitrogen_fixation', 'Nitrogen fixation of cyanos', default=.true.)
   call self%get_parameter(wdz,     'wdz',   'm/d',  'Detritus      sinking velocity', default=-4.5_rk)
   call self%get_parameter(wpz,     'wpz',   'm/d',  'Diatoms       sinking velocity', default=-0.5_rk)
   call self%get_parameter(wfz,     'wfz',   'm/d',  'Zooplankton   sinking velocity', default=0.0_rk)
   call self%get_parameter(wbz,     'wbz',   'm/d',  'Cyanobacteria sinking velocity (positive)', default=0.1_rk)
   call self%get_parameter(wpo4,    'wpo4',  'm/d',  'P-Fe in water sinking velocity', default=-1.0_rk)
   call self%get_parameter(self%sfl_po,  'sfl_po',  'mmol p/m2/d', 'constant surface phosphate flux', default=0.0015_rk)
   call self%get_parameter(self%sfl_nn,  'sfl_nn',  'mmol n/m2/d', 'constant surface nitrate flux', default=0.083_rk)
   call self%get_parameter(self%sfl_aa,  'sfl_aa',  'mmol n/m2/d','constant surface ammonium flux', default=0.06_rk)
   call self%get_parameter(self%nb,      'nb',      '1/d', 'Phytoplankton excretion rate (pl -> aa)', default=0.01_rk)
   self%nb = self%nb/secs_per_day
   call self%get_parameter(self%deltao,  'deltao',  '1/d', 'Phytoplankton mortality rate (pl -> dd)', default=0.02_rk)
   self%deltao = self%deltao/secs_per_day
   call self%get_parameter(self%nue,     'nue',     'm3/d/mmol n', 'Zooplankton respiration rate (zz -> aa)', default=0.01_rk)
   self%nue = self%nue/secs_per_day
   call self%get_parameter(self%sigma_b, 'sigma_b', 'm3/d/mmol n', 'Zooplankton mortality   rate (zz -> dd)', default=0.03_rk)
   self%sigma_b = self%sigma_b/secs_per_day
   call self%get_parameter(self%dn,      'dn',      '1/d', 'Detritus mineralization rate (dd -> aa)', default=0.003_rk)
   self%dn = self%dn/secs_per_day
   call self%get_parameter(self%dn_sed,  'dn_sed',  '1/d', 'Sediment mineralization rate (fl -> aa)', default=0.002_rk)
   self%dn_sed = self%dn_sed/secs_per_day
   call self%get_parameter(self%rp0,     'rp0',     '1/d', 'Diatoms       uptake rate', default=1.3_rk)
   self%rp0 = self%rp0/secs_per_day
   call self%get_parameter(self%rf0,     'rf0',     '1/d', 'Flagellates   uptake rate', default=0.4_rk)
   self%rf0 = self%rf0/secs_per_day
   call self%get_parameter(self%rb0,     'rb0',     '1/d', 'Cyanobacteria uptake rate', default=0.75_rk)
   self%rb0 = self%rb0/secs_per_day
   call self%get_parameter(self%cyanotll,'cyanotll','deg C', 'Cyanobacteria lower temperature limit', default=13.5_rk)
   call self%get_parameter(self%cyanosll,'cyanosll','PSU', 'Cyanobacteria lower salinity limit', default=1.0_rk)
   call self%get_parameter(self%cyanosul,'cyanosul','PSU', 'Cyanobacteria upper salinity limit', default=10.0_rk)
   call self%get_parameter(self%flagtll, 'flagtll', '(deg C)2', 'Flagellates half-saturation temp, squared', default=100.0_rk)
   call self%get_parameter(self%alphap,  'alphap',  'mmol n/m3', 'Half-saturation const, diatoms', default=0.25_rk)
   call self%get_parameter(self%alphaf,  'alphaf',  'mmol n/m3', 'Half-saturation const, flagellates', default=0.10_rk)
   call self%get_parameter(self%alphab,  'alphab',  'mmol n/m3', 'Half-saturation const, cyanobacteria', default=0.4_rk)
   call self%get_parameter(self%Yc_diat,  'Yc_diat',  'umol C/mg Chla', 'Carbon to Chla ratio (micro-mol : mg), diatoms', default=6.25_rk)
   call self%get_parameter(self%Yc_flag,  'Yc_flag',  'umol C/mg Chla', 'Carbon to Chla ratio (micro-mol : mg), flagellates', default=6.25_rk)
   call self%get_parameter(self%Yc_cyan,  'Yc_cyan',  'umol C/mg Chla', 'Carbon to Chla ratio (micro-mol : mg), cyanobacteria', default=6.25_rk)
   call self%get_parameter(self%iv,      'iv',      '1/(mmol n/m3)3', 'Ivlev constant, quadratic', default=1.2_rk)
   call self%get_parameter(self%graz,    'graz',    '1/d', 'Zooplankton grazing rate', default=0.5_rk)
   self%graz = self%graz/secs_per_day
   call self%get_parameter(self%toptz,   'toptz',   'deg C', 'Optimal temperature for grazing', default=20._rk)
   call self%get_parameter(self%zcl1,    'zcl1',    '-', 'Zooplankton closure parameter', default= 50._rk)
   call self%get_parameter(self%p0,      'p0',      'mmol n/m3', 'Diatoms       background value', default=0.001_rk)
   call self%get_parameter(self%f0,      'f0',      'mmol n/m3', 'Flagellates   background value', default=0.001_rk)
   call self%get_parameter(self%b0,      'b0',      'mmol n/m3', 'Cyanobacteria background value', default=0.001_rk)
   call self%get_parameter(self%z0,      'z0',      'mmol n/m3', 'Zooplankton   background value', default=0.001_rk)
   call self%get_parameter(self%imin_di, 'imin_di', 'W/m2', 'minimal optimal light radiation, diatoms', default=35._rk)
   call self%get_parameter(self%imin,    'imin',    'W/m2', 'minimal optimal light radiation, others', default=50._rk)
   call self%get_parameter(self%kc,      'kc',      'm2/mmol n', 'attenuation constant for self-shading effect')
   call self%get_parameter(self%q10_rec, 'q10_rec', '-', 'sediment recycling q10 rule factor', default=0.15_rk)
   call self%get_parameter(self%ade_r0,  'ade_r0',  '1/d', 'chemoautolithotrophic denitrification rate', default=0.1_rk)
   self%ade_r0 = self%ade_r0/secs_per_day
   call self%get_parameter(self%alphaade,'alphaade','mmol n/m3', 'half-saturation constant for ade', default=1.0_rk)
   call self%get_parameter(self%q10_recs,'q10_recs','-', 'sediment recycling q10 rule factor', default=0.175_rk)
   call self%get_parameter(self%rfr,     'rfr',     '-', 'Redfield ratio P/N', default=0.0625_rk)
   call self%get_parameter(self%rfc,     'rfc',     '-', 'Redfield ratio C/N', default=6.625_rk)
   call self%get_parameter(self%sedrate, 'sedrate', 'm/d', 'Detritus sedimentation rate', default=2.25_rk)
   self%sedrate = self%sedrate/secs_per_day
   call self%get_parameter(self%erorate, 'erorate', 'm/d', 'Sediment erosion rate', default=6._rk)
   self%erorate = self%erorate/secs_per_day
   call self%get_parameter(self%sedratepo4,'sedratepo4','m/d', 'P-Fe sedimentation rate', default=0.5_rk)
   self%sedratepo4 = self%sedratepo4/secs_per_day
   call self%get_parameter(self%eroratepo4,'eroratepo4','m/d', 'P-Fe     erosion rate', default=6._rk)
   self%eroratepo4 = self%eroratepo4/secs_per_day
   call self%get_parameter(self%po4ret,  'po4ret',  '-', 'Phosphate retention rate, oxic sediments', default=0.18_rk)
   call self%get_parameter(self%pburialrate,'pburialrate','-', 'Phoshate burial rate', default=0.007_rk)
   self%pburialrate = self%pburialrate/secs_per_day
   call self%get_parameter(self%fl_burialrate,'fl_burialrate','-', 'Sediment burial rate', default=0.001_rk)
   self%fl_burialrate = self%fl_burialrate/secs_per_day
   call self%get_parameter(self%pliberationrate,'pliberationrate','-', 'Phosphate liberation rate, anoxic sediments', default=0.1_rk)
   self%pliberationrate = self%pliberationrate/secs_per_day
   call self%get_parameter(self%ipo4th,  'ipo4th',  'mmol p/m2', 'Increased phosphate burial threshold', default=100._rk)
   call self%get_parameter(self%maxsed,  'maxsed',  '-', 'Maximum active sediment concentration', default=1000._rk)
   call self%get_parameter(self%br0,     'br0',     '1/d', 'Bioresuspension rate', default=0.03_rk)
   self%br0 = self%br0/secs_per_day
   call self%get_parameter(self%fds,     'fds',     '-', 'Sediment denitrification rate, oxic', default=0.7_rk)
   call self%get_parameter(self%pvel,    'pvel',    'm/d', 'Piston velocity', default=5._rk)
   self%pvel = self%pvel/secs_per_day
   call self%get_parameter(self%newflux, 'newflux', '-', 'Oxygen flux type', default=2)
   call self%get_parameter(self%tau_crit,'tau_crit','N/m2', 'Critical shear stress', default=0.07_rk)
   
   ! Register state variables
   call self%register_state_variable(self%id_pp,'pp','mmol n/m**3', 'diatoms', minimum=0.0_rk,vertical_movement=wpz/secs_per_day,no_river_dilution=.true.)
   call self%register_state_variable(self%id_ff,'ff','mmol n/m**3', 'flagellates', minimum=0.0_rk,vertical_movement=wfz/secs_per_day,no_river_dilution=.true.)
   call self%register_state_variable(self%id_bb,'bb','mmol n/m**3', 'cyanobacteria', minimum=0.0_rk,vertical_movement=wbz/secs_per_day,no_river_dilution=.true.)
   call self%register_state_variable(self%id_zz,'zz','mmol n/m**3', 'zooplankton', minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_dd,'dd','mmol n/m**3', 'detritus', minimum=0.0_rk,vertical_movement=wdz/secs_per_day,no_river_dilution=.true.)
   call self%register_state_variable(self%id_aa,'aa','mmol n/m**3', 'ammonium', minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_nn,'nn','mmol n/m**3', 'nitrate', minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_po,'po','mmol p/m**3', 'phosphate', minimum=0.0_rk,no_river_dilution=.true.)
   if (self%env_type .eq. "fresh") then
      call self%register_state_variable(self%id_o2,'o2','mmol o2/m**3','oxygen', minimum=0.0_rk,no_river_dilution=.true.)
   else
      call self%register_state_variable(self%id_o2,'o2','mmol o2/m**3','oxygen', no_river_dilution=.true.)
   end if
   call self%register_state_variable(self%id_fl,'fl','mmol n/m**2', 'fluff', minimum=0.0_rk,maximum=self%maxsed)
   call self%register_state_variable(self%id_pb,'pb','mmol p/m**2', 'PFe_s', minimum=0.0_rk)
   call self%register_state_variable(self%id_pw,'pw','mmol p/m**3', 'PFe_w', minimum=0.0_rk,vertical_movement=wpo4/secs_per_day,no_river_dilution=.true.)
     
   ! Register link to external DIC pool.
   if (self%calc_dic) call self%register_state_dependency(self%id_dic,self%dic_variable,'mmol c/m**3','dissolved inorganic carbon')

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_NO3,    'Nit',      'mgNO3N/m**3','nitrate conc in mass unit')
   call self%register_diagnostic_variable(self%id_NH4,    'Amm',      'mgNH4N/m**3','ammonium  conc in nitrogen mass unit')
   call self%register_diagnostic_variable(self%id_TN,     'TN',       'mgTN/m**3',  'total nitrogen conc in nitrogen mass unit')
   call self%register_diagnostic_variable(self%id_PO4,    'Pho',      'mgPO4P/m**3','phosphate conc in phosphorus mass unit')
   call self%register_diagnostic_variable(self%id_TP,     'TP',       'mgTP/m**3',  'total phosphorus in phosphorus mass unit')
   call self%register_diagnostic_variable(self%id_O2_mg,  'DO_mg',    'mgO2/m**3',  'oxygen in O2 mass unit')
   call self%register_diagnostic_variable(self%id_pp_chla,'pp_chla',  'mgchla/m**3','diatoms conc in chla unit')
   call self%register_diagnostic_variable(self%id_ff_chla,'ff_chla',  'mgchla/m**3','flagellates conc in chla unit')
   call self%register_diagnostic_variable(self%id_bb_chla,'bb_chla',  'mgchla/m**3','cyanobacteria conc in chla unit')
   call self%register_diagnostic_variable(self%id_tot_chla,'tot_chla','mgchla/m**3','total chlorophyll a')
   call self%register_diagnostic_variable(self%id_dPAR,    'PAR',     'W/m**2',     'photosynthetically active radiation')
   call self%register_diagnostic_variable(self%id_GPP,     'GPP',     'mmol/m**3',  'gross primary production')
   call self%register_diagnostic_variable(self%id_NCP,     'NCP',     'mmol/m**3',  'net community production')
   call self%register_diagnostic_variable(self%id_PPR,     'PPR',     'mmol/m**3/d','gross primary production rate')
   call self%register_diagnostic_variable(self%id_NPR,     'NPR',     'mmol/m**3/d','net community production rate')
   call self%register_diagnostic_variable(self%id_NFX,     'NFX',     'mmol/m**3/d','nitrogen fixation rate')
   call self%register_diagnostic_variable(self%id_DNP,     'DNP',     'mmol/m**3/d','denitrification pelagic')
   call self%register_diagnostic_variable(self%id_DNB,     'DNB',     'mmol/m**2/d','denitrification benthic')
   call self%register_diagnostic_variable(self%id_SBR,     'SBR',     'mmol/m**2',  'sediment burial')
   call self%register_diagnostic_variable(self%id_PBR,     'PBR',     'mmol/m**2/d','phosphorus burial')
   call self%register_diagnostic_variable(self%id_OFL,     'OFL',     'mmol/m**2/d','oxygen surface flux')

   ! Register environmental dependencies
   call self%register_dependency(self%id_par,  standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
   call self%register_dependency(self%id_I_0,  standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_wind, standard_variables%wind_speed)
   call self%register_dependency(self%id_taub, standard_variables%bottom_stress)

   END subroutine initialize
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of ergom1 model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class(type_msi_ergom1), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
    real(rk)           :: pp,ff,bb,zz,aa,nn,po,dd,o2,pw,par,I_0,zz0,ntemp,ntemp_eps,ptemp,phyto
    real(rk)           :: ppg,ffg,bbg,nlim,plim,rp,rf,rb
    real(rk)           :: iopt,iopt_di,ppi,ppi_di,temp,tempq,salt,ldn,nf,lpn,lpd,lp
    real(rk)           :: part_uptake,uptake_p,uptake_f,otemp,food,food_eps,lz,graz_z,ldn_N,ldn_O,ade,lzn,lzd,gg
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
    _GET_(self%id_pp,pp) ! diatoms
    _GET_(self%id_ff,ff) ! flagellates
    _GET_(self%id_bb,bb) ! cyanobacteria
    _GET_(self%id_zz,zz) ! zooplankton
    _GET_(self%id_dd,dd) ! detritus
    _GET_(self%id_aa,aa) ! ammonium
    _GET_(self%id_nn,nn) ! nitrate
    _GET_(self%id_po,po) ! phosphate
    _GET_(self%id_o2,o2) ! oxygen
    _GET_(self%id_pw,pw) ! P-Fe water

    _GET_(self%id_par,par)
    _GET_HORIZONTAL_(self%id_I_0,I_0)
    _GET_(self%id_temp,temp)
    _GET_(self%id_salt,salt)

    ! Light acclimation formulation based on surface light intensity
    iopt = max(0.5_rk*par,self%imin)
    ppi = par/iopt*exp(1.0_rk-par/iopt)

    iopt_di = max(0.5_rk*par,self%imin_di)
    ppi_di= par/iopt_di*exp(1.0_rk-par/iopt_di)

    ! if (temp .le. 0.1_rk) then
    !    ppi=0.01_rk*ppi
    !    ppi_di=0.01_rk*ppi_di
    ! end if

    tempq = max(temp,0.0_rk)
    tempq = tempq*tempq

    zz0 = self%zcl1*zz*zz

    ntemp = nn + aa
    ntemp_eps = ntemp + epsilon
    ntemp = ntemp *ntemp
    ptemp = po * po
    phyto = pp + ff + bb
    ppg = pp + self%p0
    ffg = ff + self%f0
    bbg = bb + self%b0

    ! Diatom uptake
    nlim = ntemp / (self%alphap * self%alphap + ntemp) ! MiMe eq. for IN
    plim = ptemp / (self%alphap * self%alphap * self%rfr * self%rfr + ptemp) ! MiMe eq. for IP
    rp = min (nlim, plim, ppi_di) ! Applying Liebig's minimum law: IN,IP,light
    rp = rp * self%rp0            ! Actual uptake rate

    ! Flagellate uptake
    nlim = ntemp / (self%alphaf * self%alphaf + ntemp) ! MiMe eq. for IN
    plim = ptemp / (self%alphaf * self%alphaf * self%rfr * self%rfr + ptemp) ! MiMe eq. for IP
    rf = min(nlim, plim, ppi) ! Liebig's minimum law
    rf = rf * self%rf0 * (1.0_rk + tempq / (self%flagtll + tempq)) ! Uptake rate is temp dependent

    ! Cyanobacteria uptake
    if (self%nitrogen_fixation) then
       nlim = 1.0_rk
    else
       nlim = ntemp / (self%alphab * self%alphab + ntemp) ! MiMe eq. for IN
    end if
    plim = ptemp / (self%alphab * self%alphab * self%rfr * self%rfr + ptemp) ! MiMe for IP
    rb = min(nlim, plim, ppi) !Liebig's law: IP,light. No dependence on IN for cyanos.
    rb = rb * self%rb0 / (1.0_rk + exp(self%cyanotll - temp)) !Uptake rate temp dependent
    !rb = rb / (1.0_rk + exp (salt - self%cyanosul)) / (1.0_rk + exp (self%cyanosll - salt))
    !Uptake rate is also salinity dependent for cyanos

    lpn = self%nb     ! Phytoplankton excretion rate
    lpd = self%deltao ! Phytoplankton mortality rate

    if (o2 .le. 0.0_rk) then
      ! No growth under anoxic conditions
      rp = 0.0_rk
      rf = 0.0_rk
      rb = 0.0_rk
      ! Higher mortality and no respiration if no oxygen present
      lpd = 10._rk * lpd
      lpn = 0.0_rk
    end if

    lp = lpn + lpd ! Phytoplankton loss rate (mortality + respiration)
    ! Uptake rates * concentrations / (sum of nitrate and ammonium)
    part_uptake = (rp * ppg + rf *ffg) / ntemp_eps
    uptake_p = rp * ppg / ntemp_eps
    uptake_f = rf * ffg / ntemp_eps
  !!!! NITRIFICATION RATE !!!!
    otemp = max(0.0_rk, o2) !for oxygen dependent process if o2<0 then o2=0
    ! Nitrification rate depends on oxygen availability and temperature
    nf = otemp / (0.01_rk + otemp) * 0.1_rk * exp (0.11_rk * temp)/secs_per_day
  !!!! ZOOPLANKTON RATES !!!!
    ! Food available for zooplankton. Notice lower preference for cyanos
    food = max(pp,0.0_rk) + max(ff,0.0_rk) + 0.5_rk * max(bb,0.0_rk) + self%p0 + self%b0 + self%f0
    food_eps = max(food, epsilon) ! Be sure food is positive
    gg = self%graz * (1.0_rk - exp(self%iv * food * food * (-1.0_rk))) !Grazing rate

    lzd = self%sigma_b ! Zooplankton mortality rate
    lzn = self%nue     ! Zooplankton respiration rate

    if (o2 .le. 0.0_rk) then
      ! In anoxic conditions:
      gg = 0.0_rk                ! No grazing
      lzd = 10._rk *self%sigma_b ! Higher zooplankton mortality
      lzn = 0.0_rk               ! No respiration
    endif

    lz = lzn + lzd ! Zooplankton loss rate (mortality + respiration)
    ! Zooplankton grazing depends on food availability and temperature
    graz_z = gg * (zz + self%z0)/food_eps * (1.0_rk + 2.7183_rk/self%toptz/self%toptz * tempq * exp(1.0_rk - temp * 2.0_rk / self%toptz))

    ldn = self%dn * exp (self%q10_rec*temp) ! Mineralization rate depends on temperature

    if (o2 .le. 0.0_rk .and. nn > 0.0_rk) then
      ldn_N = 5.3_rk * nn * nn / (0.001_rk + nn * nn)               ! Denitrification rate depends on nitrate availability
      ldn_O = 6.625_rk * (1.0_rk - nn * nn / (0.001_rk + nn * nn))  ! Oxygen loss due to denitrification
      ade = self%ade_r0 * nn * nn / (self%alphaade +  nn * nn) * nn ! ade rate nitrate dependent
   else
      ldn_N = 0.0_rk
      ldn_O = 6.625_rk  ! Negative oxygen = H2S production
      ade = 0.0_rk
   endif

  _ADD_SOURCE_(self%id_o2, part_uptake * (nn * 8.625_rk + aa * 6.625_rk) + rb * bbg * 6.625_rk - ldn_O * ldn * dd - 6.625_rk * (lpn * phyto + lzn * zz0) - 2.0_rk * nf * aa + ade * 0.375_rk)
  _ADD_SOURCE_(self%id_aa, ldn * dd - (part_uptake + nf) * aa + lpn * phyto + lzn * zz0)
  _ADD_SOURCE_(self%id_nn, nf * aa - part_uptake * nn - ldn_N * ldn * dd - ade)
  _ADD_SOURCE_(self%id_po, self%rfr * ( ldn * dd - (part_uptake * (aa + nn) + rb * bbg) + lpn * phyto + lzn * zz0))
  _ADD_SOURCE_(self%id_pp, uptake_p * (aa + nn) - (graz_z + lp) * pp)
  _ADD_SOURCE_(self%id_ff, uptake_f * (aa + nn) - (graz_z + lp) * ff)
  _ADD_SOURCE_(self%id_bb, rb * bbg - (0.5_rk * graz_z + lp) * bb)
  _ADD_SOURCE_(self%id_zz, graz_z * (pp + ff + 0.5_rk * bb) - lz * zz0)
  _ADD_SOURCE_(self%id_dd, lpd * phyto + lzd * zz0 - ldn * dd)

  if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic, self%rfc*( lpn * phyto + lzn * zz0 + ldn * dd - rp * ppg - rf * ffg - rb * bbg))

  _SET_DIAGNOSTIC_(self%id_NO3, nn * n_molar_mass)
  _SET_DIAGNOSTIC_(self%id_NH4, aa * n_molar_mass)
  _SET_DIAGNOSTIC_(self%id_TN, nn * n_molar_mass + aa * n_molar_mass + phyto * n_molar_mass)
  _SET_DIAGNOSTIC_(self%id_PO4, po * p_molar_mass)
  _SET_DIAGNOSTIC_(self%id_TP, po * p_molar_mass + phyto * self%rfr * p_molar_mass)
  _SET_DIAGNOSTIC_(self%id_O2_mg, o2 * o2_molar_mass)
  _SET_DIAGNOSTIC_(self%id_pp_chla, (pp * 9.5 + 2.1)/self%Yc_diat) ! relation between carbon and nitrogen from Hecky et al 1993. The stoichiometry of carbon, nitrogen, and phosphorus in particulate matter of lakes and oceans. Limnology and Oceanography, 38: 709-724.
  _SET_DIAGNOSTIC_(self%id_ff_chla, (ff * 9.5 + 2.1)/self%Yc_flag)
  _SET_DIAGNOSTIC_(self%id_bb_chla, (bb * 9.5 + 2.1)/self%Yc_cyan)
  _SET_DIAGNOSTIC_(self%id_tot_chla, (pp * 9.5 + 2.1)/self%Yc_diat + (ff * 9.5 + 2.1)/self%Yc_flag + (bb * 9.5 + 2.1)/self%Yc_cyan)
  _SET_DIAGNOSTIC_(self%id_dPAR,par)
  _SET_DIAGNOSTIC_(self%id_GPP ,rp * ppg + rf * ffg)
  _SET_DIAGNOSTIC_(self%id_NCP ,rp * ppg + rf * ffg + rb * bbg - lpn * phyto)
  _SET_DIAGNOSTIC_(self%id_PPR ,(rp * ppg + rf * ffg + rb * bbg) * secs_per_day)
  _SET_DIAGNOSTIC_(self%id_NPR ,(rp * ppg + rf * ffg + rb * bbg - lpn*phyto) * secs_per_day)
  _SET_DIAGNOSTIC_(self%id_NFX, rb * bbg * secs_per_day)
  _SET_DIAGNOSTIC_(self%id_DNP, (ldn_N * ldn *dd + ade) * secs_per_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

   END subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_msi_ergom1), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                     :: pp,ff,bb,dd
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_pp,pp) ! diatoms
   _GET_(self%id_ff,ff) ! flagellates
   _GET_(self%id_bb,bb) ! cyanobacteria
   _GET_(self%id_dd,dd) ! detritus

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(self%p0+self%f0+self%b0+pp+ff+bb+dd))

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
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
   class (type_msi_ergom1),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: pwb,pb,fl,ppb,aab,nnb,pob,ddb,oxb,taub,temp,biores
   real(rk)                   :: llds,llsd,bpds,bpsd,lldiat,recs,ldn_O,ldn_N,plib,oxlim
   real(rk)                   :: fracdenitsed,pret,pbr
   real(rk)                   :: tau_frac
   real(rk),parameter :: secs_per_day = 86400._rk
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
 !   if (self%fluff) then
   _GET_(self%id_aa,aab)
   _GET_(self%id_dd,ddb)
   _GET_(self%id_pp,ppb)
   _GET_(self%id_nn,nnb)
   _GET_(self%id_po,pob)
   _GET_(self%id_o2,oxb)
   _GET_(self%id_pw,pwb)
   _GET_HORIZONTAL_(self%id_fl,fl)
   _GET_HORIZONTAL_(self%id_pb,pb)

   _GET_HORIZONTAL_(self%id_taub,taub)
   _GET_(self%id_temp,temp)

   !increased phosphorus burial
   pbr = max(pb, pb * (pb -self%ipo4th + 1.0_rk)) 

   !bio-resuspension rate
   biores = self%br0 * max(0.0_rk, oxb) * max(0.0_rk, oxb) &
            / (max(0.0_rk, oxb) * max(0.0_rk, oxb) + 0.03_rk)  

   ! Resuspension-sedimentation rate are computed as in GOTM-BIO
! Diatoms is assumed to become detritus/fluff as soon as it settles to bottom sediments,
   ! and diatoms can therefore not be resuspended from benthic layer
  
   tau_frac = (taub-self%tau_crit)/self%tau_crit
   ! if actual taub is greater than critical tau, then do resuspension
   if (taub .gt. self%tau_crit) then
      llds=0.0_rk
      llsd=self%erorate*(taub-self%tau_crit)/self%tau_crit
      bpds=0.0_rk
      bpsd=self%eroratepo4*(taub-self%tau_crit)/self%tau_crit
      lldiat= 0.0_rk
   else
      llds=self%sedrate*(self%tau_crit-taub)/self%tau_crit
      llsd=0.0_rk
      bpds=self%sedratepo4*(self%tau_crit-taub)/self%tau_crit
      bpsd=0.0_rk
      lldiat=self%sedrate*(self%tau_crit-taub)/self%tau_crit
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
      fracdenitsed = self%fds     ! denitrification is sediments
      plib = 0.0_rk               ! no phosphate release
      pret = self%po4ret          ! phosphate is stored
   else
      fracdenitsed = 0.0_rk       ! no denitrification in sediments
      plib = self%pliberationrate ! phosphorus is liberated
      pret = 0.0_rk               ! no phosphate retention
   endif

   oxlim = max (0.0_rk,oxb) * max (0.0_rk,oxb) / (0.01_rk + max(0.0_rk,oxb) * max(0.0_rk,oxb))

   ! Sediment resuspension, detritus settling, diatom settling, bio-resuspension, mineralization and burial
   _ADD_BOTTOM_SOURCE_(self%id_fl,-llsd * fl + llds * ddb + lldiat * ppb - biores * fl - recs * fl - fl * self%fl_burialrate * fl/self%maxsed)
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
   ! Diatom settling
   _ADD_BOTTOM_FLUX_(self%id_pp, -lldiat * ppb)

   if (_AVAILABLE_(self%id_dic)) _ADD_BOTTOM_FLUX_(self%id_dic, self%rfc*recs * fl)

   ! BENTHIC DIAGNOSTIC VARIABLES
   _SET_BOTTOM_DIAGNOSTIC_(self%id_DNB,(ldn_N * recs * fl + fracdenitsed * recs * fl) * secs_per_day)
   _SET_BOTTOM_DIAGNOSTIC_(self%id_SBR,(fl * self%fl_burialrate * fl/self%maxsed) * secs_per_day)
   _SET_BOTTOM_DIAGNOSTIC_(self%id_PBR,(pbr * self%pburialrate * fl/self%maxsed) * secs_per_day)

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

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
!  type (type_msi_ergom), intent(in) :: self
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
! !IROUTINE: Surface fluxes for the ergom1 model
!
! !INTERFACE:
  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !INPUT PARAMETERS:
  class(type_msi_ergom1),intent(in) :: self

  _DECLARE_ARGUMENTS_DO_SURFACE_

  real(rk)           :: temp,wnd,salt,o2,nn,aa,po
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
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_HORIZONTAL_(self%id_wind,wnd)

   _GET_(self%id_o2,o2)
   _GET_(self%id_nn,nn)
   _GET_(self%id_aa,aa)
   _GET_(self%id_po,po)

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

   _ADD_SURFACE_FLUX_(self%id_nn,self%sfl_nn/secs_per_day)
   _ADD_SURFACE_FLUX_(self%id_aa,self%sfl_aa/secs_per_day)
   _ADD_SURFACE_FLUX_(self%id_po,self%sfl_po/secs_per_day)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------

  END MODULE fabm_msi_ergom1

!-----------------------------------------------------------------------

!Copyright (C) 2013 - Gennadi Lessin

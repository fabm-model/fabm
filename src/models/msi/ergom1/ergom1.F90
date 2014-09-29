#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!                  ERGOM (Baltic Sea Ecosystem Model)                        !
!                                                                            !
!  This version of ERGOM was coded by Gennadi Lessin (Marine Systems         !
!  Institute, Tallinn University of Technology).                             !
!  It updates the older version of ERGOM included in GOTM-BIO. This older    !
!  version is still provided in FABM under the name gotm/ergom for           !
!  historical/educational purposes.                                           !
!  The present version is based upon code provided by Thomas Neumann         !
!  for the project ECOSUPPORT, http://www.baltex-research.eu/ecosupport/.    !
!                                                                            !
!  A detailed description of the original model is given                     !
!  in T. Neumann, W. Fennel and C. Kremp, 2002. Experimental simulations     !
!  with an ecosystem model of the Baltic Sea: a nutrient load reduction      !
!  experiment. Global Biogeochemical Cycles, 16 (3).                         !
!  http://dx.doi.org/10.1029/2001GB001450.                                   !
!                                                                            !
!  The present updated version adds oxygen-dependent phosphorus dynamics     !
!  between sediment and water and the effect of bio-resuspension, as         !
!  described in T. Neumann and G. Schernewski, 2008. Eutrophication in the   !
!  Baltic Sea and shifts in nitrogen fixation analyzed with a 3D ecosystem   !
!  model, J. Mar. Sys., 74 (1.2), pp. 592-602. The present model version     !
!  also includes:                                                            !
!  a. Option for air-sea oxygen flux calculation: 1. "GOTM-BIO               !
!  type", saturation calculation according to Weiss (1970) and piston        !
!  velocity depends on wind speed, 2. Saturation using cubic polynom         !
!  approximation (HELCOM Combine manual) and 3. Calculation according to     !
!  Livingstone (1993). b. sedimentation/resuspension calculated in the       !
!  same way as it was done in GOTM-BIO version. c. Optionally DIC variable   !
!  can be included for linking with CO2 model. d. The following diagnostic   !
!  variables are calculated: PAR (photosynthetically active radiation),      !
!  GPP (gross primary production), NCP (net community production),           !
!  PPR (gross primary production rate), NPR (net primary production rate),   !
!  NFX (nitrogen fixation rate), DNP (pelagic denitrification rate),         !
!  DNB (benthic denitrification rate), SBR (sediment burial (not yet)),      !
!  PBR (phosphate burial rate), OFL (air-sea oxygen flux rate).              !
!  e. Benthic denitrification rate constant (not depth-dependent).           !
!  Updated 08/01/2013                                                        !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      type (type_diagnostic_variable_id)   :: id_dPAR,id_GPP,id_NCP,id_PPR,id_NPR,id_NFX,id_DNP
      type (type_horizontal_diagnostic_variable_id) :: id_DNB,id_SBR,id_PBR,id_OFL

      ! Model parameters
      real(rk) :: nb,deltao,nue,sigma_b,dn,dn_sed,sfl_po,sfl_aa,sfl_nn
      real(rk) :: rp0,rf0,rb0,cyanotll,cyanosll,cyanosul,flagtll
      real(rk) :: alphap,alphaf,alphab,iv,graz,toptz,zcl1,p0,f0,b0,z0
      real(rk) :: imin_di,imin,kc,q10_rec,ade_r0,alphaade,q10_recs
      real(rk) :: rfr,rfc,sedrate,erorate,sedratepo4,eroratepo4,po4ret
      real(rk) :: pburialrate,pliberationrate,ipo4th,maxsed,br0,fds,pvel,tau_crit
      integer  :: newflux
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
!   Here, the ergom1 namelist is read and the variables exported by the model are registered with FABM
!
! !INPUT PARAMETERS:
   class(type_msi_ergom1),intent(inout),target :: self
   integer,               intent(in)           :: configunit
!
! !LOCAL VARIABLES:
   real(rk)           :: pp_initial=0.001_rk    ! Initial concentration of diatoms
   real(rk)           :: ff_initial=0.001_rk    ! Initial concentration of flagellates
   real(rk)           :: bb_initial=0.001_rk    ! Initial concentration of cyanobacteria
   real(rk)           :: zz_initial=0.001_rk    ! Initial concentration of zooplankton
   real(rk)           :: dd_initial=100.0_rk    ! Initial concentration of detritus
   real(rk)           :: aa_initial=0.001_rk    ! Initial concentration of ammonium
   real(rk)           :: nn_initial=8.0_rk      ! Initial concentration of nitrate
   real(rk)           :: po_initial=1.0_rk      ! Initial concentration of phosphate
   real(rk)           :: o2_initial=300._rk     ! Initial concentration of oxygen
   real(rk)           :: fl_initial=300._rk     ! Initial concentration of sediment fluff
   real(rk)           :: pw_initial=0.001_rk    ! Initial concentration of P-Fe in water
   real(rk)           :: pb_initial=0.001_rk    ! Initial concentration of P-Fe in sediment
   real(rk)           :: sfl_po=0.0015_rk       ! Surface flux constant, phosphates
   real(rk)           :: sfl_aa=0.06_rk         ! Surface flux constant, ammonium
   real(rk)           :: sfl_nn=0.083_rk        ! Surface flux constant, nitrate
   real(rk)           :: nb=0.01_rk             ! Phytoplankton excretion rate (pl -> aa)
   real(rk)           :: deltao=0.02_rk         ! Phytoplankton mortality rate (pl -> dd)
   real(rk)           :: nue=0.01_rk            ! Zooplankton respiration rate (zz -> aa)
   real(rk)           :: sigma_b=0.03_rk        ! Zooplankton mortality   rate (zz -> dd)
   real(rk)           :: dn=0.003_rk            ! Detritus mineralization rate (dd -> aa)
   real(rk)           :: dn_sed=0.002_rk        ! Sediment mineralization rate (fl -> aa)
   real(rk)           :: rp0=1.3_rk             ! Diatoms       uptake rate
   real(rk)           :: rf0=0.4_rk             ! Flagellates   uptake rate
   real(rk)           :: rb0=0.75_rk            ! Cyanobacteria uptake rate
   real(rk)           :: cyanotll=13.5_rk       ! Cyanobacteria lower temperature limit
   real(rk)           :: cyanosll=1.0_rk        ! Cyanobacteria lower salinity limit
   real(rk)           :: cyanosul=10.0_rk       ! Cyanobacteria upper salinity limit
   real(rk)           :: flagtll=100.0_rk       ! Flagellates half-saturation temp, squared
   real(rk)           :: alphap=0.25_rk         ! Half-saturation const, diatoms
   real(rk)           :: alphaf=0.10_rk         ! Half-saturation const, flagellates
   real(rk)           :: alphab=0.4_rk          ! Half-saturation const, cyanobacteria
   real(rk)           :: iv=1.2_rk              ! Ivlev constant, quadratic
   real(rk)           :: graz=0.5_rk            ! Zooplankton grazing rate
   real(rk)           :: toptz=20._rk           ! Optimal temperature for grazing
   real(rk)           :: zcl1=50._rk            ! Zooplankton closure parameter
   real(rk)           :: p0=0.001_rk            ! Diatoms       background value
   real(rk)           :: f0=0.001_rk            ! Flagellates   backgorund value
   real(rk)           :: b0=0.001_rk            ! Cyanobacteria background value
   real(rk)           :: z0=0.001_rk            ! Zooplankton   background value
   real(rk)           :: imin_di=35.0_rk        ! minimal optimal light radiation, diatoms
   real(rk)           :: imin=50.0_rk           ! minimal optimal light radiation, others
   real(rk)           :: kc=0.03_rk             ! attenuation constant for self-shading effect
   real(rk)           :: q10_rec=0.15_rk        ! detritus recycling q10 rule factor
   real(rk)           :: ade_r0=0.1_rk          ! chemoautolithotrophic denitrification rate
   real(rk)           :: alphaade=1.0_rk        ! half-saturation constant for ade
   real(rk)           :: q10_recs=0.175_rk      ! sediment recycling q10 rule factor
   real(rk)           :: rfr=0.0625_rk          ! Redfield ratio P/N
   real(rk)           :: rfc=6.625_rk           ! Redfield ratio C/N
   real(rk)           :: wdz=-4.5_rk            ! Detritus      sinking velocity
   real(rk)           :: wpz=-0.5_rk            ! Diatoms       sinking velocity
   real(rk)           :: wfz=0.0_rk             ! Zooplankton   sinking velocity
   real(rk)           :: wbz=0.1_rk             ! Cyanobacteria sinking velocity (positive)
   real(rk)           :: wpo4=-1.0_rk           ! P-Fe in water sinking velocity
   real(rk)           :: sedrate=2.25_rk        ! Detritus sedimentation rate
   real(rk)           :: erorate=6.0_rk         ! Sediment erosion rate
   real(rk)           :: sedratepo4=0.5_rk      ! P-Fe     sedimentation rate
   real(rk)           :: eroratepo4=6.0_rk      ! P-Fe     erosion rate
   real(rk)           :: po4ret=0.18_rk         ! Phosphate retention rate,  oxic sediments
   real(rk)           :: pburialrate=0.007_rk   ! Phopshate burial rate
   real(rk)           :: pliberationrate=0.1_rk ! Phosphate liberation rate, anoxic sediments
   real(rk)           :: ipo4th=100._rk         ! Increased phosphate burial threshold
   real(rk)           :: maxsed=1000._rk        ! Maximum active sediment concentration
   real(rk)           :: br0=0.03_rk            ! Bioresuspension rate
   real(rk)           :: fds=0.7_rk             ! Sediment denitrification rate, oxic
   real(rk)           :: pvel=5._rk             ! Piston velocity
   integer            :: newflux=2              ! Oxygen flux type
                                                ! 1: Using Weiss(1970) for sat and wind dependent p_vel
                                                ! 2: Using polynom approximation S&T**3 for sat
                                                ! 3: Using simplified formulation for sat
   real(rk)           :: tau_crit=0.07_rk       ! Critical shear stress
   logical            :: calc_dic=.false.       ! Calculate dissolved organic carbon
   character(len=64)  :: dic_variable=''        ! Link variable to carbon  model
   real(rk),parameter :: secs_pr_day=86400._rk
   namelist /msi_ergom1/ pp_initial,ff_initial,bb_initial,zz_initial,  &
                        dd_initial,aa_initial,nn_initial,po_initial,   &
                        o2_initial,fl_initial,pw_initial,pb_initial,   &
                        sfl_po,sfl_aa,sfl_nn,nb,deltao,nue,sigma_b,dn, &
                        dn_sed,rp0,rf0,rb0,cyanotll,cyanosll,cyanosul, &
                        flagtll,alphap,alphaf,alphab,iv,graz,toptz,    &
                        zcl1,p0,f0,b0,z0,imin_di,imin,kc,q10_rec,      &
                        ade_r0,alphaade,q10_recs,rfr,rfc,wdz,wpz,wfz,  &
                        wbz,wpo4,sedrate,erorate,sedratepo4,eroratepo4,&
                        po4ret,pburialrate,pliberationrate,ipo4th,     &
                        maxsed,br0,fds,pvel,newflux,tau_crit,          &
                        calc_dic,dic_variable
!EOP
!-----------------------------------------------------------------------
!BOC
!  Read the namelist
   read(configunit,nml=msi_ergom1,err=99,end=100)

!  Store parameter values in our own derived type
!  NB! All rates must be provided in values per day,
!  and are converted here to values per second
   self%sfl_po   = sfl_po  ! Atmospheric nutrient fluxes are converted
   self%sfl_nn   = sfl_nn  ! to values per second in the subroutine!
   self%sfl_aa   = sfl_aa  !
   self%nb       = nb/secs_pr_day
   self%deltao   = deltao/secs_pr_day
   self%nue      = nue/secs_pr_day
   self%sigma_b  = sigma_b/secs_pr_day
   self%dn       = dn/secs_pr_day
   self%dn_sed   = dn_sed/secs_pr_day
   self%rp0      = rp0/secs_pr_day
   self%rf0      = rf0/secs_pr_day
   self%rb0      = rb0/secs_pr_day
   self%cyanotll = cyanotll
   self%cyanosll = cyanosll
   self%cyanosul = cyanosul
   self%flagtll  = flagtll
   self%alphap   = alphap
   self%alphaf   = alphaf
   self%alphab   = alphab
   self%iv       = iv
   self%graz     = graz/secs_pr_day
   self%toptz    = toptz
   self%zcl1     = zcl1
   self%p0       = p0
   self%f0       = f0
   self%b0       = b0
   self%z0       = z0
   self%imin_di  = imin_di
   self%imin     = imin
   self%kc       = kc
   self%q10_rec  = q10_rec
   self%ade_r0   = ade_r0/secs_pr_day
   self%alphaade = alphaade
   self%q10_recs = q10_recs
   self%rfr      = rfr
   self%rfc      = rfc
   self%sedrate  = sedrate/secs_pr_day
   self%erorate  = erorate/secs_pr_day
   self%sedratepo4 = sedratepo4/secs_pr_day
   self%eroratepo4 = eroratepo4/secs_pr_day
   self%po4ret   = po4ret
   self%pburialrate = pburialrate/secs_pr_day
   self%pliberationrate = pliberationrate/secs_pr_day
   self%ipo4th   = ipo4th
   self%maxsed   = maxsed
   self%br0      = br0/secs_pr_day
   self%fds      = fds
   self%pvel     = pvel/secs_pr_day
   self%newflux  = newflux
   self%tau_crit = tau_crit

   ! Register state variables
   call self%register_state_variable(self%id_pp,'pp','mmol n/m**3', 'diatoms',         pp_initial,minimum=0.0_rk,vertical_movement=wpz/secs_pr_day,no_river_dilution=.true.)
   call self%register_state_variable(self%id_ff,'ff','mmol n/m**3', 'flagellates',     ff_initial,minimum=0.0_rk,vertical_movement=wfz/secs_pr_day,no_river_dilution=.true.)
   call self%register_state_variable(self%id_bb,'bb','mmol n/m**3', 'cyanobacteria',   bb_initial,minimum=0.0_rk,vertical_movement=wbz/secs_pr_day,no_river_dilution=.true.)
   call self%register_state_variable(self%id_zz,'zz','mmol n/m**3', 'zooplankton',     zz_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_dd,'dd','mmol n/m**3', 'detritus',        dd_initial,minimum=0.0_rk,vertical_movement=wdz/secs_pr_day,no_river_dilution=.true.)
   call self%register_state_variable(self%id_aa,'aa','mmol n/m**3', 'ammonium',        aa_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_nn,'nn','mmol n/m**3', 'nitrate',         nn_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_po,'po','mmol p/m**3', 'phosphate',       po_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_o2,'o2','mmol o2/m**3','oxygen',          o2_initial,no_river_dilution=.true.)
   call self%register_state_variable(self%id_fl,'fl','mmol n/m**2', 'fluff',           fl_initial,minimum=0.0_rk,maximum=self%maxsed)
   call self%register_state_variable(self%id_pb,'pb','mmol p/m**2', 'P-Fe in sediment',pb_initial,minimum=0.0_rk)
   call self%register_state_variable(self%id_pw,'pw','mmol p/m**3', 'P-Fe in water',   pw_initial,minimum=0.0_rk,vertical_movement=wpo4/secs_pr_day,no_river_dilution=.true.)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   if (calc_dic) call self%register_state_dependency(self%id_dic,dic_variable)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2',     'photosynthetically active radiation',time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_GPP, 'GPP','mmol/m**3',  'gross primary production',           time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_NCP, 'NCP','mmol/m**3',  'net community production',           time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_PPR, 'PPR','mmol/m**3/d','gross primary production rate',      time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_NPR, 'NPR','mmol/m**3/d','net community production rate',      time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_NFX, 'NFX','mmol/m**3/d','nitrogen fixation rate',             time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_DNP, 'DNP','mmol/m**3/d','denitrification pelagic',            time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_DNB, 'DNB','mmol/m**2/d','denitrification benthic',            time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_SBR, 'SBR','mmol/m**2',  'sediment burial',                    time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_PBR, 'PBR','mmol/m**2/d','phosphorus burial',                  time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_OFL, 'OFL','mmol/m**2/d','oxygen surface flux',                time_treatment=time_treatment_averaged)

   ! Register environmental dependencies
   call self%register_dependency(self%id_par,  standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
   call self%register_dependency(self%id_I_0,  standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_wind, standard_variables%wind_speed)
   call self%register_dependency(self%id_taub, standard_variables%bottom_stress)

   return

 99 call self%fatal_error('initialize','Error reading namelist msi_ergom1.')

100 call self%fatal_error('initialize','Namelist msi_ergom1 was not found.')

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
    real(rk),parameter :: secs_pr_day=86400._rk
    real(rk),parameter :: epsilon = 0.00000000001_rk
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
    plim = ptemp / (self%alphab * self%alphab * self%rfr * self%rfr + ptemp) ! MiMe for IP
    rb = min(plim, ppi) !Liebig's law: IP,light. No dependence on IN for cyanos.
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
    nf = otemp / (0.01_rk + otemp) * 0.1_rk * exp (0.11_rk * temp)/secs_pr_day
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

  _SET_ODE_(self%id_o2, part_uptake * (nn * 8.625_rk + aa * 6.625_rk) + rb * bbg * 6.625_rk - ldn_O * ldn * dd - 6.625_rk * (lpn * phyto + lzn * zz0) - 2.0_rk * nf * aa + ade * 0.375_rk)
  _SET_ODE_(self%id_aa, ldn * dd - (part_uptake + nf) * aa + lpn * phyto + lzn * zz0)
  _SET_ODE_(self%id_nn, nf * aa - part_uptake * nn - ldn_N * ldn * dd - ade)
  _SET_ODE_(self%id_po, self%rfr * ( ldn * dd - (part_uptake * (aa + nn) + rb * bbg) + lpn * phyto + lzn * zz0))
  _SET_ODE_(self%id_pp, uptake_p * (aa + nn) - (graz_z + lp) * pp)
  _SET_ODE_(self%id_ff, uptake_f * (aa + nn) - (graz_z + lp) * ff)
  _SET_ODE_(self%id_bb, rb * bbg - (0.5_rk * graz_z + lp) * bb)
  _SET_ODE_(self%id_zz, graz_z * (pp + ff + 0.5_rk * bb) - lz * zz0)
  _SET_ODE_(self%id_dd, lpd * phyto + lzd * zz0 - ldn * dd)

  if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic, self%rfc*( lpn * phyto + lzn * zz0 + ldn * dd - rp * ppg - rf * ffg - rb * bbg))

   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP ,rp + rf)
   _SET_DIAGNOSTIC_(self%id_NCP ,rp + rf + rb - lpn * phyto)
   _SET_DIAGNOSTIC_(self%id_PPR ,(rp + rf + rb) * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR ,(rp + rf + rb - lpn*phyto) * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NFX, rb * bbg * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_DNP, (ldn_N * ldn *dd + ade) * secs_pr_day)

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
   real(rk)                   :: pwb,pb,fl,aab,nnb,pob,ddb,oxb,taub,temp,biores
   real(rk)                   :: llds,llsd,bpds,bpsd,recs,ldn_O,ldn_N,plib,oxlim
   real(rk)                   :: fracdenitsed,pret,pbr
   real(rk),parameter         :: secs_pr_day=86400.
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
   _GET_(self%id_nn,nnb)
   _GET_(self%id_po,pob)
   _GET_(self%id_o2,oxb)
   _GET_(self%id_pw,pwb)
   _GET_HORIZONTAL_(self%id_fl,fl)
   _GET_HORIZONTAL_(self%id_pb,pb)

   _GET_HORIZONTAL_(self%id_taub,taub)
   _GET_(self%id_temp,temp)

   pbr = max(pb, pb * (pb -self%ipo4th + 1.0_rk)) !increased phosphorus burial

   biores = self%br0 * max(0.0_rk, oxb) * max(0.0_rk, oxb) &
            / (max(0.0_rk, oxb) * max(0.0_rk, oxb) + 0.03_rk)  !bio-resuspension rate

   ! Resuspension-sedimentation rate are computed as in GOTM-BIO
   if (self%tau_crit .gt. taub) then
      llds=self%sedrate*(self%tau_crit-taub)/self%tau_crit
      bpds=self%sedratepo4*(self%tau_crit-taub)/self%tau_crit
   else
      llds= 0.0_rk
      bpds= 0.0_rk
   end if
   if (self%tau_crit .lt. taub) then
      llsd=self%erorate*(taub-self%tau_crit)/self%tau_crit
      bpsd=self%eroratepo4*(taub-self%tau_crit)/self%tau_crit
   else
      llsd=0.0_rk
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
      fracdenitsed = self%fds     ! denitrification is sediments
      plib = 0.0_rk               ! no phosphate release
      pret = self%po4ret          ! phosphate is stored
   else
      fracdenitsed = 0.0_rk       ! no denitrification in sediments
      plib = self%pliberationrate ! phosphorus is liberated
      pret = 0.0_rk               ! no phosphate retention
   endif

   oxlim = max (0.0_rk,oxb) * max (0.0_rk,oxb) / (0.01_rk + max(0.0_rk,oxb) * max(0.0_rk,oxb))

   ! Sediment resuspension, detritus settling, bio-resuspension and mineralization
   _SET_ODE_BEN_(self%id_fl,-llsd * fl + llds * ddb - biores * fl - recs * fl)
   ! P-Fe resuspension, sedimentation, bio-resuspension, liberation, retention and burial
   _SET_ODE_BEN_(self%id_pb,-bpsd * pb + bpds * pwb -biores * pb -plib * pb + self%rfr * recs * fl * pret * oxlim - pbr * self%pburialrate * fl/self%maxsed)

   ! Denitrification in sediments
   _SET_BOTTOM_EXCHANGE_(self%id_nn,-ldn_N * recs * fl)
   ! Oxygen consumption due to mineralization and denitrification
   _SET_BOTTOM_EXCHANGE_(self%id_o2,-ldn_O * recs * fl - 2.0_rk * fracdenitsed * recs * fl)
   ! Ammonium production due to mineralization (oxic & anoxic)
   _SET_BOTTOM_EXCHANGE_(self%id_aa,(1.0_rk - fracdenitsed) * recs * fl)
   ! Phosphate production due to mineralization (retention if oxic) and release in anoxic
   _SET_BOTTOM_EXCHANGE_(self%id_po, (1.0_rk - pret * oxlim) * self%rfr * recs * fl + plib * pb)
   ! Sediment resuspension, detritus settling, bio-resuspension
   _SET_BOTTOM_EXCHANGE_(self%id_dd, llsd * fl -llds *ddb + biores * fl)
   ! P-Fe resuspension, settling and bio-resuspension
   _SET_BOTTOM_EXCHANGE_(self%id_pw, bpsd * pb -bpds * pwb +biores * pb)

   if (_AVAILABLE_(self%id_dic)) _SET_BOTTOM_EXCHANGE_(self%id_dic, self%rfc*recs * fl)

   ! BENTHIC DIAGNOSTIC VARIABLES
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_DNB,(ldn_N * recs * fl + fracdenitsed * recs * fl) * secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_PBR,(pbr * self%pburialrate * fl/self%maxsed) * secs_pr_day)

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
  real(rk),parameter :: secs_pr_day=86400._rk
!
! !LOCAL VARIABLES:
  real(rk)                 :: p_vel,sc,flo2,osat
!  integer,parameter        :: newflux=2

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
      p_vel = p_vel/secs_pr_day
      flo2 =p_vel*(osat_weiss(temp,salt)-o2)
      _SET_SURFACE_EXCHANGE_(self%id_o2,flo2)
! Newflux=2 use polynom approximation
! to order S & T ** 3 for oxygen saturation derived from
! www.helcom.fi/Monas/CombineManual2/PartB/AnnexB-8Appendix3.pdf
   elseif (self%newflux .eq. 2) then
      osat = (10.18e0_rk + ((5.306e-3_rk - 4.8725e-5_rk * temp) *temp - 0.2785e0_rk) * temp &
          + salt * ((2.2258e-3_rk + (4.39e-7_rk * temp - 4.645e-5_rk) * temp) * temp - 6.33e-2_rk)) &
          * 44.66e0_rk
      flo2 = self%pvel * (osat - o2)
      _SET_SURFACE_EXCHANGE_(self%id_o2,flo2)
   else
      flo2 = self%pvel * (31.25_rk * (14.603_rk - 0.40215_rk * temp) - o2)
      _SET_SURFACE_EXCHANGE_(self%id_o2,flo2)
   end if

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_OFL,flo2 * secs_pr_day)

   _SET_SURFACE_EXCHANGE_(self%id_nn,self%sfl_nn/secs_pr_day)
   _SET_SURFACE_EXCHANGE_(self%id_aa,self%sfl_aa/secs_pr_day)
   _SET_SURFACE_EXCHANGE_(self%id_po,self%sfl_po/secs_pr_day)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------

  END MODULE fabm_msi_ergom1

!-----------------------------------------------------------------------

!Copyright (C) 2013 - Gennadi Lessin

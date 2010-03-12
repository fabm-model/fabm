!$Id$
#include"cppdefs.h"
#define IRON_QUOTA_LIMITATION

!-----------------------------------------------------------------------
!BOP1
!
! !MODULE: bio_npzd_fe --- simple IRON bio model \label{sec:bio_npzd_fe}
!
! !INTERFACE:
   module bio_npzd_fe
!
! !DESCRIPTION:
!  The FeNPZD model (Weber et al. 2007) which is described here consists 
!  of a simple biological NPZD-type (Oschlies and Schartau 2005) and a
!  complex iron model (Weber et al. 2005). An iron quota for
!  phytoplankton, zooplankton and detritus is introduced allowing a 
!  decoupling between the cycling of Fe and N (Weber et al. 2007). Growth 
!  of phytoplankton can depend on the iron quota by setting the define 
!  statement iron_quota_limitation.
! 
!  The iron model consists of 8 compartments (5 iron species (Fe3+, Fe2+,
!  colloidal Fe, Fe-particles, Fe-ligands), superoxide, hydrogen peroxide 
!  as well as ligands and particles. Sources and sinks of the iron species 
!  are oxidation and reduction processes (photoreduction) as well as 
!  mechanical related processes (e.g. scavenging).
!  It should be considered that the parameters for the iron model are not 
!  well known and changes might be necessary for future simulations.
!
! !USES:
!  default: all is private.
   use bio_var
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_npzd_fe, init_var_npzd_fe,                   &
          surface_fluxes_npzd_fe,light_npzd_fe, do_bio_npzd_fe, &
          clean_bio_npzd_fe
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Weber et al. 2007
!
!  $Log: bio_npzd_fe.F90,v $
!  Revision 1.2  2008-07-08 09:58:38  lars
!  adapted to changed BIO initialization algorithm
!
!  Revision 1.1  2008-03-26 08:51:44  kb
!  new directory based bio structure
!
!  Revision 1.1  2008-02-20 11:29:59  kb
!  added NPZD iron model - Weber et. all + Inga Hense
!
!
! !LOCAL VARIABLES:
!  from a namelist
   REALTYPE                  :: N_initial=4.5
   REALTYPE                  :: P_initial=0.
   REALTYPE                  :: Z_initial=0.
   REALTYPE                  :: D_initial=4.5
   REALTYPE                  :: sfl_n=0.09
   REALTYPE                  :: sfl_pt=2.38e-09
   REALTYPE                  :: sfl_f3=1.4889
   REALTYPE                  :: sfl_ho=0.0037
   REALTYPE, public          :: P0=0.0225
   REALTYPE                  :: Z0=0.0225
   REALTYPE                  :: I_min=25.
   REALTYPE                  :: Iv=1.1
   REALTYPE                  :: alpha=0.3
   REALTYPE                  :: rkw
   REALTYPE, public          :: rkc
   REALTYPE                  :: aa=0.62
   REALTYPE                  :: k_n=0.7
   REALTYPE                  :: mue_m=0.27
   REALTYPE                  :: mp = 0.04
   REALTYPE                  :: mpq= 0.025
   REALTYPE                  :: slopf= 0.925
   REALTYPE                  :: gmax= 1.575
   REALTYPE                  :: pcr= 1.6
   REALTYPE                  :: mzq= 0.34 
   REALTYPE                  :: exf= 0.01
   REALTYPE                  :: rem= 0.048
   REALTYPE                  :: w_d= -18.
   integer                   :: out_unit
   integer, parameter        :: nu=1,ph=2,zo=3,de=4

!  ---- Iron ----

   REALTYPE                  :: f3_initial=0.1
   REALTYPE                  :: f2_initial=0.1
   REALTYPE                  :: fc_initial=0.25
   REALTYPE                  :: fl_initial=0.2
   REALTYPE                  :: fp_initial=0.
   REALTYPE                  :: li_initial=0.4
   REALTYPE                  :: pt_initial=0.
   REALTYPE                  :: ho_initial=40.
   REALTYPE                  :: om_initial=0.0
   REALTYPE                  :: pf_initial=0.00033
   REALTYPE                  :: zf_initial=0.000033
   REALTYPE                  :: df_initial=0.00000033

   REALTYPE                  :: fe_d    = 1.32
   REALTYPE                  :: fe_ro2m=12960.0
   REALTYPE                  :: fe_col= 20.16
   REALTYPE                  :: fe_o2= 184.8
   REALTYPE                  :: fe_oo2m= 864.0
   REALTYPE                  :: fe_h2o2= 6.24
   REALTYPE                  :: fe_diss= 0.24
   REALTYPE                  :: fe_dismt= 2.64
   REALTYPE                  :: fe_lf= 172.8
   REALTYPE                  :: fe_ligand=20.0
   REALTYPE                  :: fe_c= 8.8
   REALTYPE                  :: fe_cd=0.2
   REALTYPE                  :: fe_pd=0.2
   REALTYPE                  :: fe_colag= 1224000.
   REALTYPE                  :: fe_star= 25000.0
   REALTYPE                  :: fe_n2cp= 0.000000159
   REALTYPE                  :: fe_phdet= 20.16
   REALTYPE                  :: fe_solub= 0.01
   REALTYPE                  :: fe_rh2o2= 15000.0
   REALTYPE                  :: fe2n= 0.033
   REALTYPE                  :: fe2n_min= 0.0033
   REALTYPE                  :: fe_irmax= 1978.0
   REALTYPE                  :: ta_f= 1.066
   REALTYPE                  :: fe_lrm= 86.4
   REALTYPE                  :: fe_so2m= 1036.8
   REALTYPE                  :: fe_cuox= 812200
   REALTYPE                  :: fe_cured= 1382
   REALTYPE                  :: fe_cutot= 1.0
   REALTYPE                  :: w_fp= -5.787037e-05
   REALTYPE                  :: w_pt= -5.787037e-05
   REALTYPE                  :: smallm = 1.e-15
   REALTYPE                  :: k_fe = 0.2
   REALTYPE                  :: fe_ld,fe_cui,fe_cuii
  
   integer, parameter        :: f3=5,f2=6,fc=7,fl=8,fp=9
   integer, parameter        :: li=10,pt=11,ho=12,om=13
   integer, parameter        :: pf=14,zf=15,df=16

!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio_npzd_fe(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio_npzd_fe.nml} is read and memory is
!  allocated - and various variables are initialised.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: namlst
   character(len=*), intent(in)   :: fname
   integer,          intent(in)   :: unit
!
! !REVISION HISTORY:
!  Original author(s): Weber et al. 2007
!
! !LOCAL VARIABLES:
   namelist /bio_npzd_fe_nml/ numc,                                     &
                N_initial,P_initial,Z_initial,D_initial,                &
                f3_initial,f2_initial,fc_initial,fl_initial,fp_initial, &
                li_initial,pt_initial,ho_initial,om_initial,            &
                pf_initial,zf_initial,df_initial,                       &
                sfl_pt,sfl_f3,sfl_ho,surface_flux_method,               &
                P0,Z0,I_min,Iv,alpha,rkw,rkc,aa,                        &
                k_n,mue_m,mp,mpq,slopf,gmax,pcr,mzq,exf,rem,w_d,        &
                fe_d,fe_ro2m,fe_col,fe_o2,fe_oo2m,fe_h2o2,              &
                fe_diss,fe_dismt,fe_lf,fe_ligand,fe_c,fe_cd,fe_pd,      &
                fe_colag,fe_star,fe_n2cp,fe_phdet,fe_solub,             &
                fe_rh2o2,fe2n,fe_irmax,ta_f,fe_lrm,fe_so2m,             &
                fe_cuox,fe_cured,fe_cutot,w_fp,w_pt,smallm,             & 
                k_fe, fe2n_min

!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_npzd_fe'

   numc=16

   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_npzd_fe_nml,err=99)
   close(namlst)

   LEVEL3 'namelist "', fname, '" read'

   n_surface_fluxes=3

!  Conversion from day to second

   mue_m= mue_m /secs_pr_day
   mp   = mp   /secs_pr_day
   mpq  = mpq  /secs_pr_day
   gmax = gmax /secs_pr_day
   pcr  = pcr  /secs_pr_day
   mzq  = mzq  /secs_pr_day
   exf  = exf  /secs_pr_day
   rem  = rem  /secs_pr_day
   w_d  = w_d  /secs_pr_day
   w_fp = w_fp /secs_pr_day
   w_pt = w_pt /secs_pr_day
   alpha = alpha/secs_pr_day

!  Fe ----

   fe_d= fe_d     /secs_pr_day
   fe_ro2m= fe_ro2m  /secs_pr_day
   fe_col= fe_col   /secs_pr_day
   fe_o2= fe_o2    /secs_pr_day
   fe_oo2m= fe_oo2m  /secs_pr_day
   fe_h2o2= fe_h2o2  /secs_pr_day
   fe_diss= fe_diss  /secs_pr_day
   fe_dismt= fe_dismt /secs_pr_day
   fe_lf= fe_lf    /secs_pr_day
   fe_c= fe_c     /secs_pr_day
   fe_cd= fe_cd     /secs_pr_day
   fe_pd= fe_pd /secs_pr_day
   fe_colag= fe_colag /secs_pr_day
   fe_star= fe_star  /secs_pr_day
   fe_phdet= fe_phdet /secs_pr_day
   fe_lrm= fe_lrm   /secs_pr_day
   fe_so2m= fe_so2m  /secs_pr_day
   fe_cuox= fe_cuox  /secs_pr_day
   fe_cured= fe_cured /secs_pr_day
!  ---

   fe_ld   = fe_lf / fe_ligand
   fe_cui  = fe_cutot * fe_cured/(fe_cured+fe_cuox)
   fe_cuii = fe_cutot * fe_cuox/(fe_cured+fe_cuox)

!   mpfq = mpq/fe2n  
!   mzfq = mzq/fe2n
!   pcrf = pcr/fe2n
! ---

!  initialize variable descriptions

   call bio_alloc_info


   var_names(1) = 'nut'
   var_units(1) = 'mmol/m**3'
   var_long(1)  = 'dissolved inorganic nitrogen'

   var_names(2) = 'phy'
   var_units(2) = 'mmol/m**3'
   var_long(2)  = 'phytoplankton'

   var_names(3) = 'zoo'
   var_units(3) = 'mmol/m**3'
   var_long(3)  = 'zooplankton'

   var_names(4) = 'det'
   var_units(4) = 'mmol/m**3'
   var_long(4)  = 'detritus'

!  --- Fe ---
   var_names(5) = 'fe3'
   var_units(5) = 'nM'
   var_long(5)  = 'Fe(III)'

   var_names(6) = 'fe2'
   var_units(6) = 'nM'
   var_long(6)  = 'Fe(II)'

   var_names(7) = 'fec'
   var_units(7) = 'nM'
   var_long(7)  = 'colloidal Fe'

   var_names(8) = 'fel'
   var_units(8) = 'nM'
   var_long(8)  = 'organically complexed Fe'

   var_names(9) = 'fep'
   var_units(9) = 'nM'
   var_long(9)  = 'adsorbed Fe'

   var_names(10) = 'lig'
   var_units(10) = 'nM'
   var_long(10)  = 'free ligands'

   var_names(11) = 'part'
   var_units(11) = 'kg/l'
   var_long(11)  = 'inorganic particles'

   var_names(12) = 'h2o2'
   var_units(12) = 'nM'
   var_long(12)  = 'hydrogen peroxide'

   var_names(13) = 'o2m'
   var_units(13) = 'nM'
   var_long(13)  = 'superoxide'

   var_names(14) = 'fphy'
   var_units(14) = 'nM'
   var_long(14)  = 'phytoplankton fe'

   var_names(15) = 'fzoo'
   var_units(15) = 'nM'
   var_long(15)  = 'zooplankton fe'

   var_names(16) = 'fdet'
   var_units(16) = 'nM'
   var_long(16)  = 'detritus fe'


   out_unit=unit

   LEVEL3 'module initialized'
#ifdef IRON_QUOTA_LIMITATION
   LEVEL3 '  (using Droop Fe quota limitation)'
#endif

   return

98 LEVEL2 'I could not open bio_npzd_fe.nml'
   LEVEL2 'If thats not what you want you have to supply bio_npzd_fe.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_npzd_fe.nml'
   return
99 FATAL 'I could not read bio_npzd_fe.nml'
   stop 'init_bio_npzd_fe'
   end subroutine init_bio_npzd_fe
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_npzd_fe
!
! !DESCRIPTION:
!  Here, the cc and ws varibles are filled with initial conditions
!
! !USES:
   IMPLICIT NONE

! !REVISION HISTORY:
!  Original author(s): Weber et al. 2007 

! !LOCAL VARIABLES:
  integer                    :: i
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev
      cc(nu,i)=n_initial
      cc(ph,i)=p_initial
      cc(zo,i)=z_initial
      cc(de,i)=d_initial

      cc(f3,i)=f3_initial
      cc(f2,i)=f2_initial
      cc(fc,i)=fc_initial
      cc(fl,i)=fl_initial
      cc(fp,i)=fp_initial
      cc(li,i)=li_initial
      cc(pt,i)=pt_initial
      cc(ho,i)=ho_initial
      cc(om,i)=om_initial
      cc(pf,i)=pf_initial 
      cc(zf,i)=zf_initial 
      cc(df,i)=df_initial 
   end do

   do i=0,nlev
      ws(nu,i) = _ZERO_
      ws(ph,i) = _ZERO_
      ws(zo,i) = _ZERO_
      ws(de,i) = w_d

      ws(f3,i)= _ZERO_
      ws(f2,i)= _ZERO_
      ws(fc,i)= _ZERO_
      ws(fl,i)= _ZERO_
      ws(fp,i)= w_fp 
      ws(li,i)= _ZERO_
      ws(pt,i)= w_pt
      ws(ho,i)= _ZERO_
      ws(om,i)= _ZERO_
      ws(pf,i) = _ZERO_
      ws(zf,i) = _ZERO_
      ws(df,i) = w_d
   end do

   sfl = _ZERO_

   posconc(nu) = 1
   posconc(ph) = 1
   posconc(zo) = 1
   posconc(de) = 1
   posconc(f3) = 1
   posconc(f2) = 1
   posconc(fc) = 1
   posconc(fl) = 1
   posconc(fp) = 1
   posconc(li) = 1
   posconc(pt) = 1
   posconc(ho) = 1
   posconc(om) = 1
   posconc(pf) = 1
   posconc(zf) = 1
   posconc(df) = 1


!  NOTE: Positive fluxes into the sea surface must have negative sign !
   select case (surface_flux_method)
      case (-1)! absolutely nothing
      case (0) ! constant

         sfl(f3)=(sfl_f3 /secs_pr_day) * fe_solub
         sfl(pt)= sfl_pt /secs_pr_day
         sfl(ho)=(sfl_ho /secs_pr_day) * fe_rh2o2

      case (2) ! from file via sfl_read

      case default
   end select

   LEVEL3 'IRON variables initialised ...'

   return

   end subroutine init_var_npzd_fe
!EOC


!----------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the IRON model
!
! !INTERFACE
   subroutine surface_fluxes_npzd_fe(nlev)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
  integer                              :: nlev
!
! !REVISION HISTORY:
!  Original author(s): Weber et al. 2007 
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

!  NOTE: Positive fluxes into the sea surface must have negative sign !
   select case (surface_flux_method)
      case (-1)! absolutely nothing
      case (0) ! constant
         sfl(f3)=(sfl_f3 /secs_pr_day) * fe_solub
         sfl(pt)= sfl_pt /secs_pr_day * 0.001 ! unit conversion
         sfl(ho)=(sfl_ho /secs_pr_day) * fe_rh2o2
      case (2) ! from file via sfl_read
         sfl(ho) =   1.*(sfl_read(1)/secs_pr_day) * fe_rh2o2
         sfl(pt) =   1.*(sfl_read(2)/secs_pr_day) * 0.001 ! unit conversion
         sfl(f3) =   1.*(sfl_read(3)/secs_pr_day) * fe_solub
      case (3) ! sfl array filled externally - for 3D models
      case default
   end select

   return
   end subroutine surface_fluxes_npzd_fe
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light properties for the FeNPZD model
!
! !INTERFACE:
   subroutine light_npzd_fe(nlev,bioshade_feedback)
!
! !DESCRIPTION:
!  Here, the photosynthetically available radiation is calculated
!  by simply assuming that the short wave part of the total
!  radiation is available for photosynthesis. The user should make
!  sure that this is consistent with the light class given in the
!  {\tt extinct} namelist of the {\tt obs.nml} file.
!  The self-shading effect is also calculated in this subroutine,
!  which may be used to consider the effect of bio-turbidity also
!  in the temperature equation (if {\tt bioshade\_feedback} is set
!  to true in {\tt bio.nml}).
!  For details, see section \ref{sec:do-bio}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                  :: nlev
   logical, intent(in)                  :: bioshade_feedback
!
! !REVISION HISTORY:
!  Original author(s): Weber et al. 2007
!
! !LOCAL VARIABLES:
   integer                   :: i,ii
   REALTYPE                  :: zz,add

!EOP
!-----------------------------------------------------------------------
!BOC
   zz = _ZERO_
   add = _ZERO_
   do i=nlev,1,-1
      add=add+0.5*h(i)*(cc(de,i)+cc(ph,i)+P0)
      zz=zz+0.5*h(i)
      par(i)=rad(nlev)*(1.-aa)*exp(-zz/rkw)*exp(-rkc*add)
      add=add+0.5*h(i)*(cc(de,i)+cc(ph,i)+P0)
      zz=zz+0.5*h(i)
      if (bioshade_feedback) bioshade_(i)=exp(-rkc*add)
   end do
   return
   end subroutine light_npzd_fe
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of geobiochemical model
!
! !INTERFACE
   subroutine do_bio_npzd_fe(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION
!
! !USES
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   logical, INTENT(in)                 :: first
   integer, INTENT(in)                 :: numc,nlev
   REALTYPE, INTENT(in)                :: cc(1:numc,0:nlev)
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)                :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(out)                :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): Weber et al. 2007 
!
! !LOCAL VARIABLES:
   integer                    :: i,j,ci
   REALTYPE, save             :: tt=0.,dt=200. 
   REALTYPE                   :: fe_fe3red,fe_fecred,fe_fe2ox,fe_lr,fe_dr
   REALTYPE                   :: feir,fe_sca,ta
   REALTYPE                   :: mue_max,mue_par,mue_nut,mue,graz
   REALTYPE                   :: fe2n_p,fe2n_z,fe_up,mue_fe,mue_f1,mue_f2
   REALTYPE                   :: exlig,det_m

!EOP
!-----------------------------------------------------------------------
!BOC
   pp = _ZERO_
   dd = _ZERO_

   do ci=1,nlev ! loop over vertical layers

!     NPZD
!     Here the rates for growth, grazing, excretion, mortality and 
!     remineralization are calculated

!     temperature dependence
      ta      = ta_f**t(ci)
!     temperature limitation
      mue_max =  mue_m * ta
!     light limitation
      mue_par = (mue_max * alpha * par(ci)) /   &
                (mue_max**2 + (alpha*par(ci))**2)**0.5 
!     nutrient limitation
      mue_nut = (cc(nu,ci)/(k_n+cc(nu,ci))) * mue_max
      mue     =  min(mue_par,mue_nut)
#ifdef IRON_QUOTA_LIMITATION
      fe2n_p = cc(pf,ci)/cc(ph,ci)    ! actual Fe:N ratio in phytoplankton
      mue_fe  =  (fe2n_p - fe2n_min)/fe2n_p * mue_max
      mue     =  min(mue,mue_fe)
#endif

      graz = (gmax*pcr*cc(ph,ci)**2)/ (gmax + pcr*cc(ph,ci)**2)  !grazing

!     Sink terms for non-negative compartments, which appear exactly
!     as or proportional to source terms for other compartments in the NPZD 
!     model part:
!     The sloppy feeding/assimilation efficiency term is assumed to be a 
!     loss term for zooplankton grazing and a source for detritus

      dd(nu,ph,ci)= mue      * cc(ph,ci)
      dd(ph,zo,ci)= graz     * cc(zo,ci)
      dd(ph,de,ci)= mpq      * cc(ph,ci)**2 
      dd(ph,nu,ci)= mp * ta  * cc(ph,ci)
      dd(zo,nu,ci)= exf* ta  * cc(zo,ci)
      dd(zo,de,ci)= mzq      * cc(zo,ci)**2 + (1-slopf)*graz*cc(zo,ci)
      dd(de,nu,ci)= rem * ta * cc(de,ci)

!     NPZD Fe internal quota
!     Here the iron uptake and content of the phytoplankton, zooplankton 
!     and detritus is calculated.

!     = fe2n (fixed Fe:N-ratio) 
      fe2n_p = cc(pf,ci)/cc(ph,ci)
!     = fe2n (fixed Fe:N-ratio) 
      fe2n_z = cc(zf,ci)/cc(zo,ci)
      fe_up  = cc(fl,ci)
      mue_fe = mue_max * (fe_up / (fe_up + k_fe))
      mue_f1 = mue_fe * cc(pf,ci)
      mue_f2 = fe2n * mue * cc(ph,ci)

!     Sink terms for non-negative compartments, which appear exactly
!     as or proportional to source terms for other compartments in the 
!     Fe-NPZD model part:
#ifdef IRON_QUOTA_LIMITATION
      dd(fl,pf,ci)= min(mue_f1,mue_f2) 
#else  
!     (fixed Fe:N-ratio) 
      dd(fl,pf,ci)= fe2n * mue * cc(ph,ci)
#endif
      dd(pf,zf,ci)= fe2n_p * graz * cc(zo,ci)
      dd(pf,df,ci)= fe2n_p * mpq  * cc(ph,ci)**2 
      dd(pf,fl,ci)= mp  * ta      * cc(pf,ci)
      dd(zf,fl,ci)= exf * ta      * cc(zf,ci)
      dd(zf,df,ci)= fe2n_z * mzq  * cc(zo,ci)**2 + fe2n_p *  &
                    (1-slopf)*graz*cc(zo,ci)
      dd(df,fl,ci)= rem * ta      * cc(df,ci)

!     Fe-module
!     Here, the processes governing the different iron forms are calculated.
!     excess (i.e. free) ligand
      exlig = cc(li,ci)
!     unit conversion detritus (micromol/l -> g/l biomass)
      det_m = fe_n2cp * cc(de,ci)
!     unit conversion light [mycroE m^-3 s^-1]
      feir  = par(ci) * 4.6
!     photoreduc.,O2-reduc.
      fe_fe3red = fe_d*feir/fe_irmax + fe_ro2m*cc(om,ci)
!     colloid. photoreduc.
      fe_fecred = fe_col*feir/fe_irmax
!     oxid. by O2,H2O2,O2-
      fe_fe2ox  = fe_o2 + fe_h2o2*cc(ho,ci) + fe_oo2m*cc(om,ci)
      fe_lr     = fe_lrm*feir/fe_irmax
      fe_dr     = fe_phdet*feir/fe_irmax
!     scavenging 
      fe_sca    = fe_star * (det_m+cc(pt,ci))

!     Sink terms for non-negative compartments, which appear exactly
!     as or proportional to source terms for other compartments in the 
!     Fe-model part:
!     fe_scavenging
      dd(f3,fp,ci)=fe_sca*cc(f3,ci)
!     fe_fe3 reduction
      dd(f3,f2,ci)=fe_fe3red*cc(f3,ci)
!     fe_colloid formation
      dd(f3,fc,ci)=fe_c*cc(f3,ci)
!     fe_cd
      dd(fc,f3,ci)=fe_cd*cc(fc,ci)
!     fe_ligand formation
      dd(f3,fl,ci)= fe_lf*exlig*cc(f3,ci)
!     fe_fe2 oxidation
      dd(f2,f3,ci)=fe_fe2ox*cc(f2,ci)
!     fe_colloid photoreduction
      dd(fc,f2,ci)=fe_fecred*cc(fc,ci)
!     fe_ligand breakup
      dd(fl,f3,ci)=fe_ld*cc(fl,ci)
!     fe_ligand photoreduction
      dd(fl,f2,ci)=fe_lr*cc(fl,ci)
!     fe_particle photoreduction
      dd(fp,f2,ci)=fe_dr*cc(fp,ci)
!     fe_pd
      dd(fp,f3,ci)=fe_pd*cc(fp,ci)
!     fe_colloidal aggragation
      dd(fc,fp,ci)=fe_colag*cc(fc,ci)*(det_m+cc(pt,ci))
      dd(om,ho,ci)= fe_oo2m*cc(f2,ci)*cc(om,ci) + & 
                 fe_cuox*fe_cui*cc(om,ci)

!     Source terms which are exactly sinks terms of other compartments or
!     proportional to them:
      do i=1,numc
         do j=1,numc
            pp(i,j,ci)=dd(j,i,ci)
         end do
      end do

!     Sink terms for positive compartments, which do not appear 
!     as source terms for other compartments:
      dd(om,om,ci)= 2.0*fe_dismt*cc(om,ci)**2 + &
                 fe_cured*fe_cuii*cc(om,ci)+ &
                 fe_ro2m*cc(f3,ci)*cc(om,ci)
      dd(ho,ho,ci)= fe_h2o2*cc(f2,ci)*cc(ho,ci) +    &
                 fe_diss*cc(ho,ci)
      dd(li,li,ci)= fe_lf*exlig*cc(f3,ci)

!     Non-conservative source terms or source and sink terms which are 
!     stoichiometrically related to other source terms:
      pp(li,li,ci)= (fe_ld + fe_lr)*cc(fl,ci)
      pp(om,om,ci)= fe_o2*cc(f2,ci)+fe_so2m * feir/fe_irmax
      pp(ho,ho,ci)= fe_dismt*cc(om,ci)**2
      if (ci .eq. 1) then     
!        bottom nudging
         pp(fc,fc,ci)= (fc_initial-cc(fc,ci))/86400.
!        bottom nudging
         pp(fl,fl,ci)= (fl_initial-cc(fl,ci))/86400.
!        bottom nudging
         pp(nu,nu,ci)= (10.623-cc(nu,ci))/86400.

      end if

   end do ! loop over vertical layers

   return

 end subroutine do_bio_npzd_fe
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine clean_bio_npzd_fe
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Weber et al. 2007
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine clean_bio_npzd_fe
!EOC

!-----------------------------------------------------------------------

   end module bio_npzd_fe

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

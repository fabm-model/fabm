#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_metu_mnemiopsis --- Mnemiopsis (comb jelly species) model
! by Baris Salihoglu (IMS METU), adapted for FABM by Jorn Bruggeman
!
! !INTERFACE:
   module fabm_metu_mnemiopsis
!
! !DESCRIPTION:
! Black Sea comb jelly (Mnemiopsis) model based on code by Baris Salihoglu.
! This model must be coupled to a lower trophic level model that provides
! prey to the Mnemiopsis population.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman, Baris Salihoglu
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_metu_mnemiopsis
!     Variable identifiers
      type (type_state_variable_id) :: id_egb, id_jb, id_ja, id_tb, id_ta, id_adb, id_ada
      type (type_state_variable_id) :: id_food, id_foodmic, id_resptarget, id_morttarget
      type (type_dependency_id)     :: id_temp

!     Model parameters
      real(rk) :: food_scale
   contains
      procedure :: initialize
      procedure :: do
   end type
!
! !PRIVATE PARAMETERS:
   real(rk),parameter :: min_bm = 1.d-15
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio\_jellyfish.nml} is read and
!  various variables (rates and settling velocities)
!  are transformed into SI units.
!
! !INPUT PARAMETERS:
   class (type_metu_mnemiopsis),intent(inout),target :: self
   integer,                     intent(in)           :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)            :: egb_initial=0.0_rk,jb_initial=0.0_rk,ja_initial=0.0_rk, &
                          tb_initial=0.0_rk,ta_initial=0.0_rk,adb_initial=0.0_rk,ada_initial=0.0_rk, &
                          food_scale=1.0_rk
   character(len=64)   :: food_source_variable,foodmic_source_variable, &
                          respiration_target_variable, mortality_target_variable
   real(rk), parameter :: secs_pr_day = 86400.
   namelist /metu_mnemiopsis/ egb_initial,jb_initial,ja_initial, &
                         tb_initial,ta_initial,adb_initial,ada_initial, &
                         food_source_variable,foodmic_source_variable, &
                         respiration_target_variable, mortality_target_variable, &
                         food_scale
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   if (configunit>0) read(configunit,nml=metu_mnemiopsis,err=99,end=100)

   ! Register state variables
   call self%register_state_variable(self%id_egb,'egb','mg C/m**3','egg biomass',        &
                                    egb_initial,minimum=0.0_rk)
   call self%register_state_variable(self%id_jb,'jb','mg C/m**3','nauplii biomass',     &
                                    jb_initial,minimum=0.0_rk)
   call self%register_state_variable(self%id_ja,'ja','#/m**3','nauplii abundance',      &
                                    ja_initial,minimum=0.0_rk)
   call self%register_state_variable(self%id_tb,'tb','mg C/m**3','transitional biomass',&
                                    tb_initial,minimum=0.0_rk)
   call self%register_state_variable(self%id_ta,'ta','#/m**3','transitional abundance', &
                                    ta_initial,minimum=0.0_rk)
   call self%register_state_variable(self%id_adb,'adb','mg C/m**3','adult biomass',      &
                                    adb_initial,minimum=0.0_rk)
   call self%register_state_variable(self%id_ada,'ada','#/m**3','adult abundance',       &
                                    ada_initial,minimum=0.0_rk)

   ! Register external state variable dependencies
   call self%register_state_dependency(self%id_food,food_source_variable)
   call self%register_state_dependency(self%id_foodmic,foodmic_source_variable)
   call self%register_state_dependency(self%id_resptarget,respiration_target_variable)
   call self%register_state_dependency(self%id_morttarget,mortality_target_variable)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)

   self%food_scale = food_scale

   return

99 call self%fatal_error('metu_mnemiopsis_init','Error reading namelist metu_mnemiopsis')
100 call self%fatal_error('metu_mnemiopsis_init','Namelist metu_mnemiopsis was not found')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Mnemiopsis model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_metu_mnemiopsis), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman, Baris Salihoglu
!
! !LOCAL VARIABLES:
   ! Internal states
   real(rk) :: egb_mn,jb_mn,ja_mn,tb_mn,ta_mn,adb_mn,ada_mn

   ! External food variables
   real(rk) :: food,foodmic,foodno

   ! Environment
   real(rk) :: temp

   ! Rates
   real(rk) :: tab_mn,teb_mn,tea_mn,tjb_mn,tja_mn,ttb_mn,tta_mn
   real(rk) :: gj_mn,mj_mn,lj_mn,gt_mn,mt_mn,lt_mn,ga_mn,ma_mn,la_mn,mea_mn

   ! Temporary variables
   real(rk) :: resp_mn
   real(rk) :: eppley
   real(rk) :: mm,mr,ma,tr,p7
   real(rk) :: AEj,AEt,AEa

   ! Model parameters
   real(rk),parameter :: betaj_mn = 0.4
   real(rk),parameter :: betat_mn = 0.09
   real(rk),parameter :: betaa_mn = 0.08
   real(rk),parameter :: mmj_mn   = 0.15 ! molting mass (max) these are after kremer 1976
   real(rk),parameter :: mmt_mn   = 1.5
   real(rk),parameter :: mma_mn   = 3.10
   real(rk),parameter :: mne_mn   = 0.0001 ! egg mg C mean mass table I
   real(rk),parameter :: mrj_mn   = 0.13 ! mg C reference  mass table I when thisincreases juvenile biomass and abundance increases
   real(rk),parameter :: mrt_mn   = 1.2
   real(rk),parameter :: mra_mn   = 2.8
   real(rk),parameter :: a        = 0.063 !Q10 factor used in the egg hatching
   real(rk), parameter :: secs_pr_day = 86400.

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Obtain current values for environmental variables.
   _GET_(self%id_temp,temp)

   ! Obtain current values for state variables
   _GET_(self%id_egb,egb_mn)
   _GET_(self%id_jb, jb_mn)
   _GET_(self%id_ja, ja_mn)
   _GET_(self%id_tb, tb_mn)
   _GET_(self%id_ta, ta_mn)
   _GET_(self%id_adb,adb_mn)
   _GET_(self%id_ada,ada_mn)
   !write (*,*) egb_mn,jb_mn,ja_mn,tb_mn,ta_mn,adb_mn,ada_mn

   ! Obtain current prey densities.
   ! These are provided externally, i.e., by lower trophic level model or data file.
   _GET_(self%id_food   ,food)
   _GET_(self%id_foodmic,foodmic)

   ! Scale external prey densities to internal unit (mg C/m**3)
   food    = food   *self%food_scale
   foodmic = foodmic*self%food_scale
   foodno  = max(1.0_rk,food/0.0024) !0.0024 is mg C per copepod

   !write (*,*) _LOCATION_,numc,food,foodmic
   !write (*,*) egb_mn,jb_mn,ja_mn,tb_mn,ta_mn,adb_mn,ada_mn

   resp_mn = 0.04+0.11*(foodno/(30.+foodno))
   lj_mn   = resp_mn+resp_mn/2.
   lt_mn   = resp_mn+resp_mn/2.
   la_mn   = resp_mn+resp_mn/2.

   eppley = exp(0.05*temp)

   ! Calculate transfer rate from egg to juveniles
   if (temp.lt.4.0) then
      teb_mn = 0.
      tea_mn = 0.
   else
      teb_mn = 0.27*exp(a*(temp-4.)) !biomass 0.27 is set so that at 25 degrees eggs hatch in 1 day
      tea_mn = teb_mn/mne_mn
   end if

   ! Calculate egg mortality
   mea_mn = 0.98+0.00002*ja_mn**2.
   if (egb_mn.lt.0.0001) mea_mn = 0.

   ! Calculate nauplii mortality
   mj_mn = 0.93+0.0003*ja_mn**2.

   if (jb_mn.gt.0.0.and.ja_mn.gt.0.0) then
      ma = jb_mn/ja_mn ! Biomass per individual nauplius
      mm = mmj_mn
      mr = mrj_mn

      if (ma.ge.mr) then
         p7=0.5
         call trans(mm,mr,ma,tr,p7)
      else
         tr=0.
      end if
      AEj=0.70!0.85-0.09*alog(foodno)
      gj_mn=betaj_mn*(((jb_mn/0.574/ja_mn)**0.574)*12.3)+0.1 !grazing of microzoo

      if (ja_mn.lt.0.5.and.ma.ge.mm) then
         gj_mn=0.
      end if

                             !              gj_mn=gj_mn*AEj*73.*food(it)!x73 try 58 kremer 1976is to convert carbon weight to dry weight
      gj_mn = gj_mn*AEj*foodmic*ja_mn/jb_mn ! ja_mn/jb_mn is to estimate ind/mgC
      gj_mn = min(gj_mn,4.d0)!after Sorokin
      tjb_mn = gj_mn*tr     !transfer rate from nauplii to copepods
      tja_mn = tjb_mn/mm    !
   else                   !if biomass or abundance is zero
      gj_mn  = 0.
      tjb_mn = 0.
      tja_mn = 0.
      lj_mn  = 0.
      mj_mn  = 0.
      jb_mn  = 0.0
      ja_mn  = 0.0

      ! Jorn: cannot influence state variables in GOTM except through temporal derivative
      !x(2)   = 0.0 !biomass
      !x(3)   = 0.00000001 !0000001!number
   end if

   ! Calculate transfer rate from "trans" to adult
   mt_mn = 0.3+0.003*ta_mn**2.

   if (ta_mn.le.0.0) then
      ma = 0.0
   else
      ma = tb_mn/ta_mn
   end if

   AEt = 0.85-0.09*log(foodno)

   if (tb_mn.gt.0.0.and.ta_mn.gt.0.0.and.Ma.le.1.0) then

      mm = mmt_mn
      mr = mrt_mn

      if (ma.ge.mr) then
         p7 = 4.
         call trans(mm,mr,ma,tr,p7)
      else
         tr = 0.
      endif

      gt_mn = betat_mn*eppley*((tb_mn/0.574/ta_mn)**(-0.5))*(0.01*foodno)**(0.65*tb_mn/ta_mn-0.65)
      ! cl_t=gt_mn
      gt_mn = gt_mn*AEt*73.*food

      ttb_mn = gt_mn*tr       ! t to a
      tta_mn = ttb_mn/mm
   end if

  !Estimate trans stage gt than 1.0
   if (tb_mn.gt.0.0.and.ta_mn.gt.0.0.and.Ma.gt.1.0) then

      Ma = tb_mn/ta_mn
      mm = mmt_mn
      mr = mrt_mn

      if (ma.ge.mr) then
         p7 = 4.
         call trans(mm,mr,ma,tr,p7)
      else
         tr = 0.
      endif

      gt_mn = betat_mn*eppley*((tb_mn/0.574/ta_mn)**(-0.5))
      gt_mn = gt_mn*AEt*73.*food
      ! to protect ma going to infinity at low food concentrations
      if (ta_mn.lt.0.5.and.ma.ge.mm) gt_mn=0.

      ttb_mn = gt_mn*tr       !c1 to c2

      tta_mn = ttb_mn/mm
   end if

   if (tb_mn.le.0.0.or.ta_mn.le.0.0) then
      gt_mn  = 0.
      ttb_mn = 0.
      tta_mn = 0.
      lt_mn  = 0.
      mt_mn  = 0.
      ma     = 0.0
      tr     = 0.0
      ! f =0.
      ! tr=0.
   end if

   ! Estimating adults
   if (adb_mn .gt. 0.0) then
      ma_mn = .01
   else
      ma_mn  = 0.0
      la_mn  = 0.0
   end if

   if (adb_mn.gt.0.0.and.ada_mn.gt.0.0) then
      Ma=adb_mn/ada_mn
      mm=mma_mn
      mr=mra_mn

      if (ma.ge.mr) then
         p7=4.
         call trans(mm,mr,ma,tr,p7)
      else
         tr=0.
      end if

      AEa = 0.85-0.09*log(foodno)
      ga_mn = betaa_mn*eppley*((adb_mn/0.574/ada_mn)**(-0.5))
      ! cl_a=ga_mn
      ga_mn = ga_mn*AEa*73.*food

      ! Calculate transfer from adults to eggs
      if (temp.lt.14.0) then
         tab_mn = 0.
      else
         tab_mn = 0.01*exp(0.115*adb_mn/0.574/ada_mn) & !*ga_mn!*tr
                  & *(0.125*25-1.875)
         ! tab_mn = 0.0   ! no reproduction case
      end if
   else !if biomass or abundance is zero
      ga_mn  = 0.
      tab_mn = 0.
      la_mn  = 0.
      ma_mn  = 0.
      ma     = 0.
   end if

   ! Increment/decrement derivatives for our own state variables
   _ADD_SOURCE_(self%id_egb,(tab_mn*adb_mn -       mea_mn       *egb_mn - teb_mn*egb_mn)/secs_pr_day) ! egg biomass
   _ADD_SOURCE_(self%id_jb ,(teb_mn*egb_mn + (gj_mn-mj_mn-lj_mn)* jb_mn - tjb_mn* jb_mn)/secs_pr_day) ! nauplii biomass
   _ADD_SOURCE_(self%id_ja ,(tea_mn*egb_mn -        mj_mn       * ja_mn - tja_mn* jb_mn)/secs_pr_day) ! nauplii abundance
   _ADD_SOURCE_(self%id_tb ,(tjb_mn* jb_mn + (gt_mn-mt_mn-lt_mn)* tb_mn - ttb_mn* tb_mn)/secs_pr_day) ! c1 biomass
   _ADD_SOURCE_(self%id_ta ,(tja_mn* jb_mn -        mt_mn       * ta_mn - tta_mn* tb_mn)/secs_pr_day) ! c1 abundance
   _ADD_SOURCE_(self%id_adb,(ttb_mn* tb_mn + (ga_mn-ma_mn-la_mn)*adb_mn - tab_mn*adb_mn)/secs_pr_day) ! adult biomass
   _ADD_SOURCE_(self%id_ada,(tta_mn* tb_mn         -ma_mn       *ada_mn                )/secs_pr_day) ! adult abundance

   ! Increment/decrement derivatives for external variables (prey, respiration/mortality target variables)
   _ADD_SOURCE_(self%id_foodmic,   -(gj_mn*jb_mn)/self%food_scale/secs_pr_day)
   _ADD_SOURCE_(self%id_food,      -(gt_mn*tb_mn+ga_mn*adb_mn)/self%food_scale/secs_pr_day)
   _ADD_SOURCE_(self%id_resptarget, (lj_mn*jb_mn + lt_mn*tb_mn + la_mn*adb_mn)/self%food_scale/secs_pr_day)
   _ADD_SOURCE_(self%id_morttarget, (mea_mn*egb_mn + mj_mn*jb_mn + mt_mn*tb_mn + ma_mn*adb_mn)/self%food_scale/secs_pr_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

   pure subroutine trans(mm,mr,ma,t,p7)
      ! here mm is max and ma actual and mr is the reference weight
      real(rk),intent(in ) :: mm,mr,ma,p7
      real(rk),intent(out) :: t

      ! p7 =4.
      t = ((ma-mr)**p7)/((ma-mr)**p7+(mm-mr)**p7)
   end subroutine trans


!-----------------------------------------------------------------------

   end module fabm_metu_mnemiopsis

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

!$Id$
#include "rmbm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: rmbm_mnemiopsis --- Mnemiopsis (comb jelly species) model
! by Baris Salihoglu (IMS METU), adapted for RMBM by Jorn Bruggeman
!
! !INTERFACE:
   module rmbm_mnemiopsis
!
! !DESCRIPTION:
! Black Sea comb jelly (Mnemiopsis) model based on code by Baris Salihoglu.
! This model must be coupled to a lower trophic level model that provides
! prey to the Mnemiopsis population.
!
! !USES:
   use rmbm_types
   use rmbm_driver

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_mnemiopsis, mnemiopsis_init, mnemiopsis_do
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman, Baris Salihoglu
!
!
! !PUBLIC DERIVED TYPES:
   type type_mnemiopsis
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_ :: id_egb, id_jb, id_ja, id_tb, id_ta, id_adb, id_ada
      _TYPE_STATE_VARIABLE_ID_ :: id_food, id_foodmic, id_resptarget, id_morttarget
      _TYPE_DEPENDENCY_ID_     :: id_temp
      
!     Model parameters
      REALTYPE :: food_scale
   end type
!
! !PRIVATE PARAMETERS:
   REALTYPE,parameter :: min_bm = 1.d-15
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine mnemiopsis_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the bio namelist {\tt bio\_jellyfish.nml} is read and 
!  various variables (rates and settling velocities) 
!  are transformed into SI units.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_mnemiopsis),intent(out)   :: self
   type (type_model_info),intent(inout) :: modelinfo
   integer,               intent(in )   :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE            :: egb_initial=_ZERO_,jb_initial=_ZERO_,ja_initial=_ZERO_, &
                          tb_initial=_ZERO_,ta_initial=_ZERO_,adb_initial=_ZERO_,ada_initial=_ZERO_, &
                          food_scale=_ONE_
   character(len=64)   :: food_source_variable,foodmic_source_variable, &
                          respiration_target_variable, mortality_target_variable
   REALTYPE, parameter :: secs_pr_day = 86400.
   namelist /mnemiopsis/ egb_initial,jb_initial,ja_initial, &
                         tb_initial,ta_initial,adb_initial,ada_initial, &
                         food_source_variable,foodmic_source_variable, &
                         respiration_target_variable, mortality_target_variable, &
                         food_scale
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=mnemiopsis,err=99)

   ! Register state variables
   self%id_egb = register_state_variable(modelinfo,'egb','mg C/m**3','egg biomass',        &
                                    egb_initial,minimum=_ZERO_)
   self%id_jb  = register_state_variable(modelinfo,'jb','mg C/m**3','nauplii biomass',     &
                                    jb_initial,minimum=_ZERO_)
   self%id_ja  = register_state_variable(modelinfo,'ja','#/m**3','nauplii abundance',      &
                                    ja_initial,minimum=_ZERO_)
   self%id_tb  = register_state_variable(modelinfo,'tb','mg C/m**3','transitional biomass',&
                                    tb_initial,minimum=_ZERO_)
   self%id_ta  = register_state_variable(modelinfo,'ta','#/m**3','transitional abundance', &
                                    ta_initial,minimum=_ZERO_)
   self%id_adb = register_state_variable(modelinfo,'adb','mg C/m**3','adult biomass',      &
                                    adb_initial,minimum=_ZERO_)
   self%id_ada = register_state_variable(modelinfo,'ada','#/m**3','adult abundance',       &
                                    ada_initial,minimum=_ZERO_)
                                    
   ! Register external state variable dependencies
   self%id_food       = register_state_dependency(modelinfo,food_source_variable)
   self%id_foodmic    = register_state_dependency(modelinfo,foodmic_source_variable)
   self%id_resptarget = register_state_dependency(modelinfo,respiration_target_variable)
   self%id_morttarget = register_state_dependency(modelinfo,mortality_target_variable)
   
   ! Register environmental dependencies
   self%id_temp = register_dependency(modelinfo, varname_temp)

   self%food_scale = food_scale

   return

99 call fatal_error('mnemiopsis_init','Error reading namelist mnemiopsis')
   
   end subroutine mnemiopsis_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Mnemiopsis model
!
! !INTERFACE:
   _PURE_ subroutine mnemiopsis_do(self,_RMBM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_mnemiopsis), intent(in) :: self
   _DECLARE_RMBM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman, Baris Salihoglu
!
! !LOCAL VARIABLES:
   ! Internal states
   REALTYPE :: egb_mn,jb_mn,ja_mn,tb_mn,ta_mn,adb_mn,ada_mn
   
   ! External food variables
   REALTYPE :: food,foodmic,foodno
   
   ! Environment
   REALTYPE :: temp
   
   ! Rates
   REALTYPE :: tab_mn,teb_mn,tea_mn,tjb_mn,tja_mn,ttb_mn,tta_mn
   REALTYPE :: gj_mn,mj_mn,lj_mn,gt_mn,mt_mn,lt_mn,ga_mn,ma_mn,la_mn,mea_mn
   
   ! Temporary variables
   REALTYPE :: resp_mn
   REALTYPE :: eppley
   REALTYPE :: mm,mr,ma,tr,p7
   REALTYPE :: AEj,AEt,AEa

   ! Model parameters
   REALTYPE,parameter :: betaj_mn = 0.4
   REALTYPE,parameter :: betat_mn = 0.09
   REALTYPE,parameter :: betaa_mn = 0.08
   REALTYPE,parameter :: mmj_mn   = 0.15 ! molting mass (max) these are after kremer 1976
   REALTYPE,parameter :: mmt_mn   = 1.5
   REALTYPE,parameter :: mma_mn   = 3.10
   REALTYPE,parameter :: mne_mn   = 0.0001 ! egg mg C mean mass table I
   REALTYPE,parameter :: mrj_mn   = 0.13 ! mg C reference  mass table I when thisincreases juvenile biomass and abundance increases
   REALTYPE,parameter :: mrt_mn   = 1.2
   REALTYPE,parameter :: mra_mn   = 2.8
   REALTYPE,parameter :: a        = 0.063 !Q10 factor used in the egg hatching
   REALTYPE, parameter :: secs_pr_day = 86400.

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _RMBM_LOOP_BEGIN_

   ! Obtain current values for environmental variables
   _GET_DEPENDENCY_(self%id_temp,temp)

   ! Obtain current values for state variables
   _GET_STATE_(self%id_egb,egb_mn)
   _GET_STATE_(self%id_jb, jb_mn)
   _GET_STATE_(self%id_ja, ja_mn)
   _GET_STATE_(self%id_tb, tb_mn)
   _GET_STATE_(self%id_ta, ta_mn)
   _GET_STATE_(self%id_adb,adb_mn)
   _GET_STATE_(self%id_ada,ada_mn)
   !write (*,*) egb_mn,jb_mn,ja_mn,tb_mn,ta_mn,adb_mn,ada_mn

   _GET_STATE_(self%id_food   ,food)
   _GET_STATE_(self%id_foodmic,foodmic)
   food = food*self%food_scale
   foodmic = foodmic*self%food_scale
   foodno  = max(_ONE_,food/0.0024) !0.0024 is mg C per copepod
   
   !write (*,*) LOCATION,numc,food,foodmic
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
   _SET_ODE_(self%id_egb,(tab_mn*adb_mn -       mea_mn       *egb_mn - teb_mn*egb_mn)/secs_pr_day) ! egg biomass
   _SET_ODE_(self%id_jb ,(teb_mn*egb_mn + (gj_mn-mj_mn-lj_mn)* jb_mn - tjb_mn* jb_mn)/secs_pr_day) ! nauplii biomass
   _SET_ODE_(self%id_ja ,(tea_mn*egb_mn -        mj_mn       * ja_mn - tja_mn* jb_mn)/secs_pr_day) ! nauplii abundance
   _SET_ODE_(self%id_tb ,(tjb_mn* jb_mn + (gt_mn-mt_mn-lt_mn)* tb_mn - ttb_mn* tb_mn)/secs_pr_day) ! c1 biomass
   _SET_ODE_(self%id_ta ,(tja_mn* jb_mn -        mt_mn       * ta_mn - tta_mn* tb_mn)/secs_pr_day) ! c1 abundance
   _SET_ODE_(self%id_adb,(ttb_mn* tb_mn + (ga_mn-ma_mn-la_mn)*adb_mn - tab_mn*adb_mn)/secs_pr_day) ! adult biomass
   _SET_ODE_(self%id_ada,(tta_mn* tb_mn         -ma_mn       *ada_mn                )/secs_pr_day) ! adult abundance
   
   ! Increment/decrement derivatives for external variables (prey, respiration/mortality target variables)
   _SET_ODE_(self%id_foodmic,   -(gj_mn*jb_mn)/self%food_scale/secs_pr_day)
   _SET_ODE_(self%id_food,      -(gt_mn*tb_mn+ga_mn*adb_mn)/self%food_scale/secs_pr_day)
   _SET_ODE_(self%id_resptarget, (lj_mn*jb_mn + lt_mn*tb_mn + la_mn*adb_mn)/self%food_scale/secs_pr_day)
   _SET_ODE_(self%id_morttarget, (mea_mn*egb_mn + mj_mn*jb_mn + mt_mn*tb_mn + ma_mn*adb_mn)/self%food_scale/secs_pr_day)

   ! Leave spatial loops (if any)
   _RMBM_LOOP_END_

   end subroutine mnemiopsis_do
!EOC

   _PURE_ subroutine trans(mm,mr,ma,t,p7)
      implicit none
      
      ! here mm is max and ma actual and mr is the reference weight 
      REALTYPE,intent(in ) :: mm,mr,ma,p7
      REALTYPE,intent(out) :: t
      
      ! p7 =4.
      t = ((ma-mr)**p7)/((ma-mr)**p7+(mm-mr)**p7)
   end subroutine trans


!-----------------------------------------------------------------------

   end module rmbm_mnemiopsis

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

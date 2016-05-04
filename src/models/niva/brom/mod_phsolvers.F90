!
!    Copyright 2013 Guy Munhoven
!
!    This file is part of SolveSAPHE.

!    SolveSAPHE is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SolveSAPHE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with SolveSAPHE.  If not, see <http://www.gnu.org/licenses/>.
!





! **********************
! Precompiler directives
! **********************

! VARIANT_BACASTOWORIG:
! - if not defined, use secant iterations on [H+]
! - if defined, use secant iterations on X = SQRT(K_1 K_2)/[H+]

#undef VARIANT_BACASTOWORIG
!#define DEBUG_PHSOLVERS

! **************************
! End precompiler directives
! **************************



MODULE MOD_PHSOLVERS

use fabm_types, only: rk

IMPLICIT NONE

! General parameters

REAL(KIND=rk), PARAMETER :: pp_rdel_ah_target = 1.E-8_rk
REAL(KIND=rk), PARAMETER :: pp_ln10 = 2.302585092994045684018_rk


! Maximum number of iterations for each method

INTEGER, PARAMETER :: jp_maxniter_atgen    = 50
INTEGER, PARAMETER :: jp_maxniter_icacfp   = 50
INTEGER, PARAMETER :: jp_maxniter_bacastow = 50
INTEGER, PARAMETER :: jp_maxniter_atsec    = 50
INTEGER, PARAMETER :: jp_maxniter_ocmip    = 50
INTEGER, PARAMETER :: jp_maxniter_atfast   = 50

! Bookkeeping variables for each method

! - SOLVE_AT_GENERAL
INTEGER :: niter_atgen    = jp_maxniter_atgen

! - SOLVE_AT_ICACFP
INTEGER :: niter_icacfp   = jp_maxniter_icacfp

! - SOLVE_AT_BACASTOW
INTEGER :: niter_bacastow = jp_maxniter_bacastow

! - SOLVE_AT_GENERAL_SEC
INTEGER :: niter_atsec    = jp_maxniter_atsec

! - SOLVE_AT_OCMIP
INTEGER :: niter_ocmip    = jp_maxniter_ocmip

! - SOLVE_AT_FAST (variant of SOLVE_AT_GENERAL w/o bracketing
INTEGER :: niter_atfast   = jp_maxniter_atfast


! Keep the following functions private to avoid conflicts with
! other modules that provide similar ones.

PRIVATE AHINI_FOR_AT, SOLVE_AC


CONTAINS

!===============================================================================
 SUBROUTINE ANW_INFSUP(p_dictot, p_bortot,                                     &
                       p_po4tot, p_siltot,  p_nh4tot, p_h2stot,                &
                       p_so4tot, p_flutot,                                     &
                       p_alknw_inf, p_alknw_sup)
!===============================================================================

! Subroutine returns the lower and upper bounds of "non-water-selfionization"
! contributions to total alkalinity (the infimum and the supremum), i.e
! inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])

IMPLICIT NONE


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)  :: p_dictot
REAL(KIND=rk), INTENT(IN)  :: p_bortot
REAL(KIND=rk), INTENT(IN)  :: p_po4tot
REAL(KIND=rk), INTENT(IN)  :: p_siltot
REAL(KIND=rk), INTENT(IN)  :: p_nh4tot
REAL(KIND=rk), INTENT(IN)  :: p_h2stot
REAL(KIND=rk), INTENT(IN)  :: p_so4tot
REAL(KIND=rk), INTENT(IN)  :: p_flutot
REAL(KIND=rk), INTENT(OUT) :: p_alknw_inf
REAL(KIND=rk), INTENT(OUT) :: p_alknw_sup


!==============================================================================


! p_alknw_inf = -\Sum_i m_i Xtot_i

! p_alknw_inf =-p_dictot*0._rk &          ! n = 2, m = 0
!              -p_bortot*0._rk &          ! n = 1, m = 0
!              -p_po4tot*1._rk &          ! n = 3, m = 1
!              -p_siltot*0._rk &          ! n = 1, m = 0
!              -p_nh4tot*0._rk &          ! n = 1, m = 0
!              -p_h2stot*0._rk &          ! n = 1, m = 0
!              -p_so4tot*1._rk &          ! n = 1, m = 1
!              -p_flutot*1._rk            ! n = 1, m = 1

p_alknw_inf =    -p_po4tot - p_so4tot - p_flutot


! p_alknw_sup = \Sum_i (n_i - m_i) Xtot_i

! p_alknw_sup = p_dictot*(2._rk-0._rk) &  ! n = 2, m = 0
!               p_bortot*(1._rk-0._rk) &  ! n = 1, m = 0
!               p_po4tot*(3._rk-1._rk) &  ! n = 3, m = 1
!               p_siltot*(1._rk-0._rk) &  ! n = 1, m = 0
!               p_nh4tot*(1._rk-0._rk) &  ! n = 1, m = 0
!               p_h2stot*(1._rk-0._rk) &  ! n = 1, m = 0
!               p_so4tot*(1._rk-1._rk) &  ! n = 1, m = 1
!               p_flutot*(1._rk-1._rk)    ! n = 1, m = 1

p_alknw_sup =   p_dictot + p_dictot + p_bortot &
              + p_po4tot + p_po4tot + p_siltot &
              + p_nh4tot + p_h2stot

RETURN

!===============================================================================
 END SUBROUTINE ANW_INFSUP
!===============================================================================




!===============================================================================
 FUNCTION EQUATION_AT(p_alktot, p_h,       p_dictot, p_bortot,                 &
                      p_po4tot, p_siltot,  p_nh4tot, p_h2stot,                 &
                      p_so4tot, p_flutot,                                      &
                      p_deriveqn)
!===============================================================================

USE MOD_CHEMCONST


IMPLICIT NONE

REAL(KIND=rk) :: EQUATION_AT


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)            :: p_alktot
REAL(KIND=rk), INTENT(IN)            :: p_h
REAL(KIND=rk), INTENT(IN)            :: p_dictot
REAL(KIND=rk), INTENT(IN)            :: p_bortot
REAL(KIND=rk), INTENT(IN)            :: p_po4tot
REAL(KIND=rk), INTENT(IN)            :: p_siltot
REAL(KIND=rk), INTENT(IN)            :: p_nh4tot
REAL(KIND=rk), INTENT(IN)            :: p_h2stot
REAL(KIND=rk), INTENT(IN)            :: p_so4tot
REAL(KIND=rk), INTENT(IN)            :: p_flutot
REAL(KIND=rk), INTENT(OUT), OPTIONAL :: p_deriveqn


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk) :: znumer_dic, zdnumer_dic, zdenom_dic, zalk_dic, zdalk_dic
REAL(KIND=rk) :: znumer_bor, zdnumer_bor, zdenom_bor, zalk_bor, zdalk_bor
REAL(KIND=rk) :: znumer_po4, zdnumer_po4, zdenom_po4, zalk_po4, zdalk_po4
REAL(KIND=rk) :: znumer_sil, zdnumer_sil, zdenom_sil, zalk_sil, zdalk_sil
REAL(KIND=rk) :: znumer_nh4, zdnumer_nh4, zdenom_nh4, zalk_nh4, zdalk_nh4
REAL(KIND=rk) :: znumer_h2s, zdnumer_h2s, zdenom_h2s, zalk_h2s, zdalk_h2s
REAL(KIND=rk) :: znumer_so4, zdnumer_so4, zdenom_so4, zalk_so4, zdalk_so4
REAL(KIND=rk) :: znumer_flu, zdnumer_flu, zdenom_flu, zalk_flu, zdalk_flu
REAL(KIND=rk) ::                                      zalk_wat, zdalk_wat


!==============================================================================


! H2CO3 - HCO3 - CO3 : n=2, m=0
znumer_dic = 2._rk*api2_dic + p_h*       api1_dic
zdenom_dic =       api2_dic + p_h*(      api1_dic + p_h)
zalk_dic   = p_dictot * (znumer_dic/zdenom_dic)

! B(OH)3 - B(OH)4 : n=1, m=0
znumer_bor =       api1_bor
zdenom_bor =       api1_bor + p_h
zalk_bor   = p_bortot * (znumer_bor/zdenom_bor)

! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
znumer_po4 = 3._rk*api3_po4 + p_h*(2._rk*api2_po4 + p_h* api1_po4)
zdenom_po4 =       api3_po4 + p_h*(      api2_po4 + p_h*(api1_po4 + p_h))
zalk_po4   = p_po4tot * (znumer_po4/zdenom_po4 - 1._rk) ! Zero level of H3PO4 = 1

! H4SiO4 - H3SiO4 : n=1, m=0
znumer_sil =       api1_sil
zdenom_sil =       api1_sil + p_h
zalk_sil   = p_siltot * (znumer_sil/zdenom_sil)

! NH4 - NH3 : n=1, m=0
znumer_nh4 =       api1_nh4
zdenom_nh4 =       api1_nh4 + p_h
zalk_nh4   = p_nh4tot * (znumer_nh4/zdenom_nh4)

! H2S - HS : n=1, m=0
znumer_h2s =       api1_h2s
zdenom_h2s =       api1_h2s + p_h
zalk_h2s   = p_h2stot * (znumer_h2s/zdenom_h2s)

! HSO4 - SO4 : n=1, m=1
znumer_so4 =       api1_so4
zdenom_so4 =       api1_so4 + p_h
zalk_so4   = p_so4tot * (znumer_so4/zdenom_so4 - 1._rk)

! HF - F : n=1, m=1
znumer_flu =       api1_flu
zdenom_flu =       api1_flu + p_h
zalk_flu   = p_flutot * (znumer_flu/zdenom_flu - 1._rk)

! H2O - OH
zalk_wat   = api1_wat/p_h - p_h/aphscale


EQUATION_AT =    zalk_dic + zalk_bor + zalk_po4 + zalk_sil &
               + zalk_nh4 + zalk_h2s + zalk_so4 + zalk_flu &
               + zalk_wat - p_alktot


IF(PRESENT(p_deriveqn)) THEN

   ! H2CO3 - HCO3 - CO3 : n=2
   zdnumer_dic = api1_dic*api2_dic + p_h*(4._rk*api2_dic                       &
                                   + p_h*       api1_dic)
   zdalk_dic   = -p_dictot*(zdnumer_dic/zdenom_dic**2)

   ! B(OH)3 - B(OH)4 : n=1
   zdnumer_bor = api1_bor
   zdalk_bor   = -p_bortot*(zdnumer_bor/zdenom_bor**2)

   ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3
   zdnumer_po4 = api2_po4*api3_po4 + p_h*(4._rk*api1_po4*api3_po4              &
                                   + p_h*(9._rk*api3_po4 + api1_po4*api2_po4   &
                                   + p_h*(4._rk*api2_po4                       &
                                   + p_h*       api1_po4)))
   zdalk_po4   = -p_po4tot * (zdnumer_po4/zdenom_po4**2)

   ! H4SiO4 - H3SiO4 : n=1
   zdnumer_sil = api1_sil
   zdalk_sil   = -p_siltot * (zdnumer_sil/zdenom_sil**2)

   ! NH4 - NH3 : n=1
   zdnumer_nh4 = api1_nh4
   zdalk_nh4   = -p_nh4tot * (zdnumer_nh4/zdenom_nh4**2)

   ! H2S - HS : n=1
   zdnumer_h2s = api1_h2s
   zdalk_h2s   = -p_h2stot * (zdnumer_h2s/zdenom_h2s**2)

   ! HSO4 - SO4 : n=1
   zdnumer_so4 = api1_so4
   zdalk_so4   = -p_so4tot * (zdnumer_so4/zdenom_so4**2)

   ! HF - F : n=1
   zdnumer_flu = api1_flu
   zdalk_flu   = -p_flutot * (zdnumer_flu/zdenom_flu**2)

   p_deriveqn =   zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
                + zdalk_nh4 + zdalk_h2s + zdalk_so4 + zdalk_flu &
                - api1_wat/p_h**2 - 1._rk/aphscale
ENDIF

RETURN

!===============================================================================
 END FUNCTION EQUATION_AT
!===============================================================================




!===============================================================================
 SUBROUTINE AHINI_FOR_AT(p_alkcb, p_dictot, p_bortot, p_hini)
!===============================================================================

! Subroutine returns the root for the 2nd order approximation of the
! DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic polynomial)
! around the local minimum, if it exists.

! Returns * 1E-03_rk if p_alkcb <= 0
!         * 1E-10_rk if p_alkcb >= 2*p_dictot + p_bortot
!         * 1E-07_rk if 0 < p_alkcb < 2*p_dictot + p_bortot
!                    and the 2nd order approximation does not have a solution


USE MOD_CHEMCONST, ONLY : api1_dic, api2_dic, api1_bor

IMPLICIT NONE


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)             ::  p_alkcb, p_dictot, p_bortot
REAL(KIND=rk), INTENT(OUT)            ::  p_hini


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk)  ::  zca, zba
REAL(KIND=rk)  ::  zd, zsqrtd, zhmin
REAL(KIND=rk)  ::  za2, za1, za0


!==============================================================================


IF (p_alkcb <= 0._rk) THEN
  p_hini = 1.e-3_rk
ELSEIF (p_alkcb >= (2._rk*p_dictot + p_bortot)) THEN
  p_hini = 1.e-10_rk
ELSE
  zca = p_dictot/p_alkcb
  zba = p_bortot/p_alkcb

  ! Coefficients of the cubic polynomial
  za2 = api1_bor*(1._rk - zba) + api1_dic*(1._rk-zca)
  za1 = api1_dic*api1_bor*(1._rk - zba - zca) + api2_dic*(1._rk - (zca+zca))
  za0 = api2_dic*api1_bor*(1._rk - zba - (zca+zca))


                                        ! Taylor expansion around the minimum

  zd = za2*za2 - 3._rk*za1              ! Discriminant of the quadratic equation
                                        ! for the minimum close to the root

  IF(zd > 0._rk) THEN                   ! If the discriminant is positive

    zsqrtd = SQRT(zd)

    IF(za2 < 0) THEN
      zhmin = (-za2 + zsqrtd)/3._rk
    ELSE
      zhmin = -za1/(za2 + zsqrtd)
    ENDIF

    p_hini = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)

  ELSE

    p_hini = 1.e-7_rk

  ENDIF

ENDIF


RETURN

!===============================================================================
 END SUBROUTINE AHINI_FOR_AT
!===============================================================================





!===============================================================================
 FUNCTION SOLVE_AT_GENERAL(p_alktot, p_dictot, p_bortot,                       &
                           p_po4tot, p_siltot, p_nh4tot, p_h2stot,             &
                           p_so4tot, p_flutot, p_hini,   p_val)
!===============================================================================

! Universal pH solver that converges from any given initial value,
! determines upper an lower bounds for the solution if required

USE MOD_CHEMCONST, ONLY: api1_wat, aphscale

IMPLICIT NONE

REAL(KIND=rk) :: SOLVE_AT_GENERAL


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)            :: p_alktot
REAL(KIND=rk), INTENT(IN)            :: p_dictot
REAL(KIND=rk), INTENT(IN)            :: p_bortot
REAL(KIND=rk), INTENT(IN)            :: p_po4tot
REAL(KIND=rk), INTENT(IN)            :: p_siltot
REAL(KIND=rk), INTENT(IN)            :: p_nh4tot
REAL(KIND=rk), INTENT(IN)            :: p_h2stot
REAL(KIND=rk), INTENT(IN)            :: p_so4tot
REAL(KIND=rk), INTENT(IN)            :: p_flutot
REAL(KIND=rk), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=rk), INTENT(OUT), OPTIONAL :: p_val


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk)  ::  zh_ini, zh, zh_prev, zh_lnfactor
REAL(KIND=rk)  ::  zalknw_inf, zalknw_sup
REAL(KIND=rk)  ::  zh_min, zh_max
REAL(KIND=rk)  ::  zdelta, zh_delta
REAL(KIND=rk)  ::  zeqn, zdeqndh, zeqn_absmin

LOGICAL        :: l_exitnow
REAL(KIND=rk), PARAMETER :: pz_exp_threshold = 1.0_rk


!==============================================================================


IF(PRESENT(p_hini)) THEN

   zh_ini = p_hini

ELSE

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_GENERAL] Calling AHINI_FOR_AT for h_ini'
#endif

   CALL AHINI_FOR_AT(p_alktot, p_dictot, p_bortot, zh_ini)

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_GENERAL] h_ini :', zh_ini
#endif


ENDIF

 CALL ANW_INFSUP(p_dictot, p_bortot,                                           &
                 p_po4tot, p_siltot,  p_nh4tot, p_h2stot,                      &
                 p_so4tot, p_flutot,                                           &
                 zalknw_inf, zalknw_sup)

zdelta = (p_alktot-zalknw_inf)**2 + 4._rk*api1_wat/aphscale

IF(p_alktot >= zalknw_inf) THEN
   zh_min = 2._rk*api1_wat /( p_alktot-zalknw_inf + SQRT(zdelta) )
ELSE
   zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._rk
ENDIF


zdelta = (p_alktot-zalknw_sup)**2 + 4._rk*api1_wat/aphscale

IF(p_alktot <= zalknw_sup) THEN
   zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._rk
ELSE
   zh_max = 2._rk*api1_wat /( p_alktot-zalknw_sup + SQRT(zdelta) )
ENDIF

#if defined(DEBUG_PHSOLVERS)
      PRINT*, '[SOLVE_AT_GENERAL] h_min :', zh_min
      PRINT*, '[SOLVE_AT_GENERAL] h_max :', zh_max
#endif

zh = MAX(MIN(zh_max, zh_ini), zh_min)
!zh = SQRT(zh_max*zh_min)              ! Uncomment this line for the
                                       ! "safe" initialisation test

niter_atgen        = 0                 ! Reset counters of iterations

zeqn_absmin        = HUGE(1._rk)


DO

   IF(niter_atgen >= jp_maxniter_atgen) THEN
      zh = -1._rk
      EXIT
   ENDIF

   zh_prev = zh

   zeqn = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot, p_nh4tot, p_h2stot,                  &
                      p_so4tot, p_flutot, P_DERIVEQN = zdeqndh)

   ! Adapt bracketing interval
   IF(zeqn > 0._rk) THEN
      zh_min = zh_prev
   ELSEIF(zeqn < 0._rk) THEN
      zh_max = zh_prev
   ELSE
      ! zh is the root; unlikely but, one never knows
      EXIT
   ENDIF


   ! Now determine the next iterate zh
   niter_atgen = niter_atgen + 1



   IF(ABS(zeqn) >= 0.5_rk*zeqn_absmin) THEN

      ! if the function evaluation at the current point is
      ! not decreasing faster than with a bisection step (at least linearly)
      ! in absolute value take one bisection step on [ph_min, ph_max]
      ! ph_new = (ph_min + ph_max)/2d0
      ! In terms of [H]_new:
      ! [H]_new = 10**(-ph_new)
      !         = 10**(-(ph_min + ph_max)/2d0)
      !         = SQRT(10**(-(ph_min + phmax)))
      !         = SQRT(zh_max * zh_min)
     
      zh = SQRT(zh_max * zh_min)

      zh_lnfactor = (zh - zh_prev)/zh_prev ! Required to test convergence below

   ELSE

      ! dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
      !           = -zdeqndh * LOG(10) * [H]
      ! \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*LOG(10))

      ! pH_new = pH_old + \deltapH

      ! [H]_new = 10**(-pH_new)
      !         = 10**(-pH_old - \Delta pH)
      !         = [H]_old * 10**(-zeqn/(zdeqndh*[H]_old*LOG(10)))
      !         = [H]_old * EXP(-LOG(10)*zeqn/(zdeqndh*[H]_old*LOG(10)))
      !         = [H]_old * EXP(-zeqn/(zdeqndh*[H]_old))

      zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

      IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
         zh          = zh_prev*EXP(zh_lnfactor)
      ELSE
         zh_delta    = zh_lnfactor*zh_prev
         zh          = zh_prev + zh_delta
      ENDIF

#if defined(DEBUG_PHSOLVERS)
      PRINT*, '[SOLVE_AT_GENERAL] testing zh :', zh, zeqn, zh_lnfactor
#endif


      IF( zh < zh_min ) THEN
         ! if [H]_new < [H]_min
         ! i.e., if ph_new > ph_max then
         ! take one bisection step on [ph_prev, ph_max]
         ! ph_new = (ph_prev + ph_max)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_max)/2d0)
         !         = SQRT(10**(-(ph_prev + phmax)))
         !         = SQRT([H]_old*10**(-ph_max))
         !         = SQRT([H]_old * zh_min)

         zh                = SQRT(zh_prev * zh_min)

         zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below

      ENDIF

      IF( zh > zh_max ) THEN
         ! if [H]_new > [H]_max
         ! i.e., if ph_new < ph_min, then
         ! take one bisection step on [ph_min, ph_prev]
         ! ph_new = (ph_prev + ph_min)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_min)/2d0)
         !         = SQRT(10**(-(ph_prev + ph_min)))
         !         = SQRT([H]_old*10**(-ph_min))
         !         = SQRT([H]_old * zhmax)

         zh                = SQRT(zh_prev * zh_max)

         zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below

      ENDIF


   ENDIF

   zeqn_absmin = MIN( ABS(zeqn), zeqn_absmin)


   ! Stop iterations once |\delta{[H]}/[H]| < rdel
   ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
   ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

   ! Alternatively:
   ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
   !             ~ 1/LOG(10) * |\Delta [H]|/[H]
   !             < 1/LOG(10) * rdel

   ! Hence |zeqn/(zdeqndh*zh)| < rdel

   ! rdel <-- pp_rdel_ah_target

   l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

   IF(l_exitnow) EXIT

ENDDO

SOLVE_AT_GENERAL = zh

IF(PRESENT(p_val)) THEN

   IF(zh > 0._rk) THEN

      p_val = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot, p_nh4tot, p_h2stot,              &
                          p_so4tot, p_flutot)

   ELSE

      p_val = HUGE(1._rk)

   ENDIF

ENDIF

RETURN

!===============================================================================
 END FUNCTION SOLVE_AT_GENERAL
!===============================================================================




!===============================================================================
 FUNCTION AC_FROM_AT(p_alktot, p_h,                 p_bortot,                  &
                     p_po4tot, p_siltot,  p_nh4tot, p_h2stot,                  &
                     p_so4tot, p_flutot)
!===============================================================================

USE MOD_CHEMCONST


IMPLICIT NONE

REAL(KIND=rk) :: AC_FROM_AT


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)  :: p_alktot
REAL(KIND=rk), INTENT(IN)  :: p_h
REAL(KIND=rk), INTENT(IN)  :: p_bortot
REAL(KIND=rk), INTENT(IN)  :: p_po4tot
REAL(KIND=rk), INTENT(IN)  :: p_siltot
REAL(KIND=rk), INTENT(IN)  :: p_nh4tot
REAL(KIND=rk), INTENT(IN)  :: p_h2stot
REAL(KIND=rk), INTENT(IN)  :: p_so4tot
REAL(KIND=rk), INTENT(IN)  :: p_flutot


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk) :: znumer_bor, zdenom_bor, zalk_bor
REAL(KIND=rk) :: znumer_po4, zdenom_po4, zalk_po4
REAL(KIND=rk) :: znumer_sil, zdenom_sil, zalk_sil
REAL(KIND=rk) :: znumer_nh4, zdenom_nh4, zalk_nh4
REAL(KIND=rk) :: znumer_h2s, zdenom_h2s, zalk_h2s
REAL(KIND=rk) :: znumer_so4, zdenom_so4, zalk_so4
REAL(KIND=rk) :: znumer_flu, zdenom_flu, zalk_flu
REAL(KIND=rk) ::                         zalk_wat


!==============================================================================


! B(OH)3 - B(OH)4 : n=1, m=0
znumer_bor =       api1_bor
zdenom_bor =       api1_bor + p_h
zalk_bor   = p_bortot * (znumer_bor/zdenom_bor)

! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
znumer_po4 = 3._rk*api3_po4 + p_h*(2._rk*api2_po4 + p_h* api1_po4)
zdenom_po4 =       api3_po4 + p_h*(      api2_po4 + p_h*(api1_po4 + p_h))
zalk_po4   = p_po4tot * (znumer_po4/zdenom_po4 - 1._rk) ! Zero level of H3PO4 = 1

! H4SiO4 - H3SiO4 : n=1, m=0
znumer_sil =       api1_sil
zdenom_sil =       api1_sil + p_h
zalk_sil   = p_siltot * (znumer_sil/zdenom_sil)

! NH4 - NH3 : n=1, m=0
znumer_nh4 =       api1_nh4
zdenom_nh4 =       api1_nh4 + p_h
zalk_nh4   = p_nh4tot * (znumer_nh4/zdenom_nh4)

! H2S - HS : n=1, m=0
znumer_h2s =       api1_h2s
zdenom_h2s =       api1_h2s + p_h
zalk_h2s   = p_h2stot * (znumer_h2s/zdenom_h2s)

! HSO4 - SO4 : n=1, m=1
znumer_so4 =       api1_so4
zdenom_so4 =       api1_so4 + p_h
zalk_so4   = p_so4tot * (znumer_so4/zdenom_so4 - 1._rk)

! HF - F : n=1, m=1
znumer_flu =       api1_flu
zdenom_flu =       api1_flu + p_h
zalk_flu   = p_flutot * (znumer_flu/zdenom_flu - 1._rk)

! H2O - OH
zalk_wat   = api1_wat/p_h - p_h/aphscale


AC_FROM_AT =   p_alktot &
             - (  zalk_bor + zalk_po4 + zalk_sil &
                + zalk_nh4 + zalk_h2s + zalk_so4 + zalk_flu &
                + zalk_wat)


RETURN

!===============================================================================
 END FUNCTION AC_FROM_AT
!===============================================================================




!===============================================================================
 FUNCTION SOLVE_AC(p_alkc, p_dictot)
!===============================================================================

! Function returns the solution of the DIC - A_C equation for [H]
! Returns -1 if A_C <= 0 or A_C >= 2*DIC

USE MOD_CHEMCONST, ONLY : api1_dic, api2_dic

IMPLICIT NONE

REAL(KIND=rk)  ::  SOLVE_AC


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)  ::  p_alkc, p_dictot


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk)  ::  zca, zsqrtdelta, za1, za0


!==============================================================================


IF((p_alkc <= 0._rk) .OR. (p_alkc >= (p_dictot+p_dictot))) THEN

  SOLVE_AC = -1._rk

ELSE

  zca = p_dictot/p_alkc
  za1 = api1_dic*(1._rk - zca)
  za0 = api2_dic*(1._rk - zca - zca)

  zsqrtdelta = SQRT(za1**2 - 4._rk*za0)

  IF(za1 > 0._rk) THEN
    SOLVE_AC = -2._rk*za0/( za1 + zsqrtdelta )
  ELSE
    SOLVE_AC = ( -za1 + zsqrtdelta )/2._rk
  ENDIF

ENDIF

RETURN

!===============================================================================
 END FUNCTION SOLVE_AC
!===============================================================================




!===============================================================================
 FUNCTION SOLVE_AT_ICACFP( p_alktot, p_dictot, p_bortot,                       &
                           p_po4tot, p_siltot, p_nh4tot, p_h2stot,             &
                           p_so4tot, p_flutot, p_hini,   p_val)
!===============================================================================

! Function returns the solution of the Alk_T-pH equation derived by
! the iterated carbonate alkalinity correction method (fixed point iteration)
! Returns -1 if the iterations did not converge (divergence or number of
! iterations exceeded.

IMPLICIT NONE

REAL(KIND=rk)  ::  SOLVE_AT_ICACFP


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)            :: p_alktot
REAL(KIND=rk), INTENT(IN)            :: p_dictot
REAL(KIND=rk), INTENT(IN)            :: p_bortot
REAL(KIND=rk), INTENT(IN)            :: p_po4tot
REAL(KIND=rk), INTENT(IN)            :: p_siltot
REAL(KIND=rk), INTENT(IN)            :: p_nh4tot
REAL(KIND=rk), INTENT(IN)            :: p_h2stot
REAL(KIND=rk), INTENT(IN)            :: p_so4tot
REAL(KIND=rk), INTENT(IN)            :: p_flutot
REAL(KIND=rk), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=rk), INTENT(OUT), OPTIONAL :: p_val


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk)  ::  zh_ini, zh, zh_prev, zalk_dic
REAL(KIND=rk)  ::  zalknw_inf, zalknw_sup
REAL(KIND=rk)  ::  zdelta
REAL(KIND=rk)  ::  zeqn, zdeqndh, zeqn_absmin

LOGICAL        ::  l_exitnow


!==============================================================================


IF(PRESENT(p_hini)) THEN

   zh_ini = p_hini

ELSE

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_ICACFP] Calling AHINI_FOR_AT for h_ini'
#endif

   CALL AHINI_FOR_AT(p_alktot, p_dictot, p_bortot, zh_ini)

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_ICACFP] h_ini :', zh_ini
#endif


ENDIF


zh = zh_ini

niter_icacfp        = 0                 ! Reset counters of iterations


DO

   niter_icacfp = niter_icacfp + 1

   IF(niter_icacfp > jp_maxniter_icacfp) THEN
      zh = -1._rk
      EXIT
   ENDIF


   zh_prev  = zh

   zalk_dic = AC_FROM_AT(p_alktot, zh,                  p_bortot,              &
                         p_po4tot, p_siltot,  p_nh4tot, p_h2stot,              &
                         p_so4tot, p_flutot)

   zh = SOLVE_AC(zalk_dic, p_dictot)

   IF(zh < 0._rk) THEN
                                       ! IF zh < 0,  the quadratic equation
                                       ! in DIC and ALK_C does not have any
                                       ! positive root.
     l_exitnow = .TRUE.

   ELSE

     l_exitnow = (ABS((zh_prev - zh)/zh) < pp_rdel_ah_target)

   ENDIF


   IF(l_exitnow) EXIT

ENDDO


SOLVE_AT_ICACFP = zh


IF(PRESENT(p_val)) THEN

   IF(zh > 0._rk) THEN

     p_val = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,               &
                         p_po4tot, p_siltot, p_nh4tot, p_h2stot,               &
                         p_so4tot, p_flutot)

   ELSE

     p_val = HUGE(1._rk)

   ENDIF

ENDIF
RETURN

!===============================================================================
END FUNCTION SOLVE_AT_ICACFP
!===============================================================================




!===============================================================================
 FUNCTION SOLVE_AT_BACASTOW(p_alktot, p_dictot, p_bortot,                      &
                            p_po4tot, p_siltot, p_nh4tot, p_h2stot,            &
                            p_so4tot, p_flutot, p_hini,   p_val)
!===============================================================================

#if defined(VARIANT_BACASTOWORIG)
USE MOD_CHEMCONST, ONLY : api2_dic
#endif

IMPLICIT NONE

! fixed point iteration

REAL(KIND=rk)  ::  SOLVE_AT_BACASTOW


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)            :: p_alktot
REAL(KIND=rk), INTENT(IN)            :: p_dictot
REAL(KIND=rk), INTENT(IN)            :: p_bortot
REAL(KIND=rk), INTENT(IN)            :: p_po4tot
REAL(KIND=rk), INTENT(IN)            :: p_siltot
REAL(KIND=rk), INTENT(IN)            :: p_nh4tot
REAL(KIND=rk), INTENT(IN)            :: p_h2stot
REAL(KIND=rk), INTENT(IN)            :: p_so4tot
REAL(KIND=rk), INTENT(IN)            :: p_flutot
REAL(KIND=rk), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=rk), INTENT(OUT), OPTIONAL :: p_val


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk)  ::  zh_ini, zh, zfh, zh_1, zfh_1, zh_2, zfh_2, zalk_dic
REAL(KIND=rk)  ::  zalknw_inf, zalknw_sup
REAL(KIND=rk)  ::  zdelta
REAL(KIND=rk)  ::  zeqn, zdeqndh, zeqn_absmin

LOGICAL        :: l_exitnow

#if defined(VARIANT_BACASTOWORIG)
REAL(KIND=rk)  ::  zscale, zx, zfx, zx_1, zfx_1, zx_2, zfx_2
#endif


!==============================================================================


IF(PRESENT(p_hini)) THEN

   zh_ini = p_hini

ELSE

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_BACASTOW] Calling AHINI_FOR_AT for h_ini'
#endif

   CALL AHINI_FOR_AT(p_alktot, p_dictot, p_bortot, zh_ini)

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_BACASTOW] h_ini :', zh_ini
#endif


ENDIF


! Prepare the secant iterations: two initial pairs
! (zh, SOLVE_AC(AC_FROM_AT(..., zh,...)) are required

! - first iterate (will become $n-2$ iterate at the first secant evaluation)

niter_bacastow     = 1                 ! Set counter of iterations

zh_1 = zh_ini
zalk_dic = AC_FROM_AT(p_alktot, zh_1,                p_bortot,                 &
                      p_po4tot, p_siltot,  p_nh4tot, p_h2stot,                 &
                      p_so4tot, p_flutot)

zfh_1 = SOLVE_AC(zalk_dic, p_dictot)

IF(zfh_1 < 0._rk) THEN
                                       ! IF zfh_1 < 0, the quadratic equation
                                       ! in DIC and ALK_C does not have any
                                       ! positive root => return to caller
   SOLVE_AT_BACASTOW = zfh_1

   IF(PRESENT(p_val)) p_val = HUGE(1._rk)

   RETURN

ENDIF


                                       ! Check if convergence criterion
                                       ! possibly already fulfilled
IF(ABS((zfh_1 - zh_1)/zh_1) < pp_rdel_ah_target) THEN

   SOLVE_AT_BACASTOW = zfh_1           ! root found (we know that zfh_1 > 0)

   IF(PRESENT(p_val)) THEN

      p_val = EQUATION_AT(p_alktot, zfh_1,    p_dictot, p_bortot,              &
                          p_po4tot, p_siltot, p_nh4tot, p_h2stot,              &
                          p_so4tot, p_flutot)

   ENDIF

   RETURN
ENDIF


! - second iterate (will become $n-1$ iterate at the first secant)

niter_bacastow     = 2                 ! Set counter of iterations

zh = zfh_1
zalk_dic = AC_FROM_AT(p_alktot, zh,                  p_bortot,                 &
                      p_po4tot, p_siltot,  p_nh4tot, p_h2stot,                 &
                      p_so4tot, p_flutot)

zfh = SOLVE_AC(zalk_dic, p_dictot)

IF(zfh < 0._rk) THEN
                                       ! IF zfh < 0, the quadratic equation
                                       ! in DIC and ALK_C does not have any
                                       ! positive root => return to caller
   SOLVE_AT_BACASTOW = zfh

   IF(PRESENT(p_val)) p_val = HUGE(1._rk)

   RETURN

ENDIF

                                       ! Check if convergence criterion
                                       ! possibly already fulfilled
IF(ABS((zfh - zh)/zh) < pp_rdel_ah_target) THEN

   SOLVE_AT_BACASTOW = zfh             ! root found (we know that zfh > 0)

   IF(PRESENT(p_val)) THEN

      p_val = EQUATION_AT(p_alktot, zfh,      p_dictot, p_bortot,              &
                          p_po4tot, p_siltot, p_nh4tot, p_h2stot,              &
                          p_so4tot, p_flutot)

   ENDIF

   RETURN

ENDIF


#if defined(VARIANT_BACASTOWORIG)
! Bacastows original method applies the secant method not on H,
! but on its scaled inverse X = SQRT(K_C1*K_C2)/H = SQRT(api2_dic)/H

zscale = SQRT(api2_dic)

zx_1  = zscale/zh_1
zfx_1 = zscale/zfh_1
zx    = zscale/zh
zfx   = zscale/zfh
#endif


DO

   niter_bacastow = niter_bacastow + 1

   IF(niter_bacastow > jp_maxniter_bacastow) THEN
      zfh = -1._rk
      EXIT
   ENDIF


#if defined(VARIANT_BACASTOWORIG)

   zx_2  = zx_1                        ! X_{n-2}
   zfx_2 = zfx_1                       ! f(X_{n-2}), and F(X_{n-2}) = X_{n-2} - f(X_{n-2})

   zx_1  = zx                          ! X_{n-1}
   zfx_1 = zfx                         ! f(X_{n-1}), and F(X_{n-1}) = X_{n-1} - f(X_{n-1})


   zx = zx_1 - ( zx_1 - zfx_1 )/(( zx_1 - zfx_1 - zx_2 + zfx_2)/( zx_1 - zx_2 ))

   zh = zscale/zx

   zalk_dic = AC_FROM_AT(p_alktot, zh,                  p_bortot,              &
                         p_po4tot, p_siltot,  p_nh4tot, p_h2stot,              &
                         p_so4tot, p_flutot)

   zfh = SOLVE_AC(zalk_dic, p_dictot)  ! evaluate f(H_{n})
   zfx = zscale/zfh                    ! scaled inverse of f(H_{n})

#else

   zh_2  = zh_1                        ! H_{n-2}
   zfh_2 = zfh_1                       ! f(H_{n-2}), and F(H_{n-2}) = H_{n-2} - f(H_{n-2})

   zh_1  = zh                          ! H_{n-1}
   zfh_1 = zfh                         ! f(H_{n-1}), and F(H_{n-1}) = H_{n-1} - f(H_{n-1})


                                       ! Calculate the iterate H_{n} by the secant method:
                                       !
                                       ! H_{n} = H_{n-1} - F(H_{n-1})*(H_{n-2} - H_{n-1})/(F(H_{n-2})-F(H_{n-1}))
                                       !
                                       !       = H_{n-1} -  (H_{n-1} - f(H_{n-1}))
                                       !                   *(H_{n-2} - H_{n-1})
                                       !                   /(H_{n-2} - f(H_{n-2}) - H_{n-1} + f(H_{n-1}))
                                       !
                                       !       = H_{n-1} -  (H_{n-1} - f(H_{n-1}))
                                       !                   *(H_{n-1} - H_{n-2})
                                       !                   /(H_{n-1} - f(H_{n-1}) - H_{n-2} + f(H_{n-2}))

   !zh = zh_1 - ( zh_1 - zfh_1 )*( zh_1 - zh_2 )/( zh_1 - zfh_1 - zh_2 + zfh_2)
   zh = zh_1 - ( zh_1 - zfh_1 )/(( zh_1 - zfh_1 - zh_2 + zfh_2)/( zh_1 - zh_2 ))

   zalk_dic = AC_FROM_AT(p_alktot, zh,                  p_bortot,              &
                         p_po4tot, p_siltot,  p_nh4tot, p_h2stot,              &
                         p_so4tot, p_flutot)

                                       ! evaluate f( H_{n})
   zfh = SOLVE_AC(zalk_dic, p_dictot)

#endif


   IF(zfh < 0._rk) THEN
                                       ! IF zfh < 0, there is no solution to the quadratic
                                       ! equation in DIC and ALK_C => return to caller
     l_exitnow = .TRUE.

   ELSE

#if defined(VARIANT_BACASTOWORIG)
     ! actually: (zfh - zh)/zh = (zx - zfx)/zfx
     l_exitnow = (ABS((zfx - zx)/zx) < pp_rdel_ah_target)
#else
     l_exitnow = (ABS((zfh - zh)/zh) < pp_rdel_ah_target)
#endif


   ENDIF


   IF(l_exitnow) EXIT

ENDDO


SOLVE_AT_BACASTOW = zfh


IF(PRESENT(p_val)) THEN

   IF(zfh > 0._rk) THEN

     p_val = EQUATION_AT(p_alktot, zfh,      p_dictot, p_bortot,               &
                         p_po4tot, p_siltot, p_nh4tot, p_h2stot,               &
                         p_so4tot, p_flutot)

   ELSE

     p_val = HUGE(1._rk)

   ENDIF

ENDIF


RETURN

!===============================================================================
END FUNCTION SOLVE_AT_BACASTOW
!===============================================================================





!===============================================================================
 FUNCTION SOLVE_AT_GENERAL_SEC(p_alktot, p_dictot, p_bortot,                   &
                               p_po4tot, p_siltot, p_nh4tot, p_h2stot,         &
                               p_so4tot, p_flutot, p_hini,   p_val)
!===============================================================================

! Universal pH solver that converges from any given initial value,
! determines upper an lower bounds for the solution if required

USE MOD_CHEMCONST, ONLY: api1_wat, aphscale

IMPLICIT NONE

REAL(KIND=rk) :: SOLVE_AT_GENERAL_SEC


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)            :: p_alktot
REAL(KIND=rk), INTENT(IN)            :: p_dictot
REAL(KIND=rk), INTENT(IN)            :: p_bortot
REAL(KIND=rk), INTENT(IN)            :: p_po4tot
REAL(KIND=rk), INTENT(IN)            :: p_siltot
REAL(KIND=rk), INTENT(IN)            :: p_nh4tot
REAL(KIND=rk), INTENT(IN)            :: p_h2stot
REAL(KIND=rk), INTENT(IN)            :: p_so4tot
REAL(KIND=rk), INTENT(IN)            :: p_flutot
REAL(KIND=rk), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=rk), INTENT(OUT), OPTIONAL :: p_val


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk)  ::  zh_ini, zh, zh_1, zh_2, zh_delta
REAL(KIND=rk)  ::  zalknw_inf, zalknw_sup
REAL(KIND=rk)  ::  zh_min, zh_max
REAL(KIND=rk)  ::  zeqn, zeqn_1, zeqn_2, zeqn_absmin
REAL(KIND=rk)  ::  zdelta

LOGICAL        ::  l_exitnow


!==============================================================================


IF(PRESENT(p_hini)) THEN

   zh_ini = p_hini

ELSE

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_GENERAL_SEC] Calling AHINI_FOR_AT for h_ini'
#endif

   CALL AHINI_FOR_AT(p_alktot, p_dictot, p_bortot, zh_ini)

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_GENERAL_SEC] h_ini :', zh_ini
#endif


ENDIF

 CALL ANW_INFSUP(p_dictot, p_bortot,                                      &
                 p_po4tot, p_siltot,  p_nh4tot, p_h2stot,                 &
                 p_so4tot, p_flutot,                                      &
                 zalknw_inf, zalknw_sup)

zdelta = (p_alktot-zalknw_inf)**2 + 4._rk*api1_wat/aphscale

IF(p_alktot >= zalknw_inf) THEN
   zh_min = 2._rk*api1_wat /( p_alktot-zalknw_inf + SQRT(zdelta) )
ELSE
   zh_min = aphscale*(-(p_alktot-zalknw_inf) + SQRT(zdelta) ) / 2._rk
ENDIF


zdelta = (p_alktot-zalknw_sup)**2 + 4._rk*api1_wat/aphscale

IF(p_alktot <= zalknw_sup) THEN
   zh_max = aphscale*(-(p_alktot-zalknw_sup) + SQRT(zdelta) ) / 2._rk
ELSE
   zh_max = 2._rk*api1_wat /( p_alktot-zalknw_sup + SQRT(zdelta) )
ENDIF

#if defined(DEBUG_PHSOLVERS)
      PRINT*, '[SOLVE_AT_GENERAL_SEC] h_min :', zh_min
      PRINT*, '[SOLVE_AT_GENERAL_SEC] h_max :', zh_max
#endif

zh = MAX(MIN(zh_max, zh_ini), zh_min)
!zh = SQRT(zh_max*zh_min)              ! Uncomment this line for the
                                       ! "safe" initialisation test

niter_atsec        = 0                 ! Reset counters of iterations




! Prepare the secant iterations: two initial (zh, zeqn) pairs are required
! We have the starting value, that needs to be completed by the evaluation
! of the equation value it produces.

! Complete the initial value with its equation evaluation
! (will take the role of the $n-2$ iterate at the first secant evaluation)

niter_atsec = 0                        ! zh_2 is the initial value;


zh_2   = zh
zeqn_2 = EQUATION_AT(p_alktot, zh_2,     p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot, p_nh4tot, p_h2stot,                 &
                     p_so4tot, p_flutot)


zeqn_absmin        = ABS(zeqn_2)


! Adapt bracketing interval and heuristically set zh_1

IF(zeqn_2 < 0._rk) THEN
                                       ! If zeqn_2 < 0, then we adjust zh_max:
                                       ! we can be sure that zh_min < zh_2 < zh_max.
   zh_max = zh_2
                                       ! for zh_1, try 25% (0.1 pH units) below the current zh_max,
                                       ! but stay above SQRT(zh_min*zh_max), which would be equivalent
                                       ! to a bisection step on [pH@zh_min, pH@zh_max]
   zh_1   = MAX(zh_max/1.25_rk, SQRT(zh_min*zh_max))

ELSEIF(zeqn_2 > 0._rk) THEN
                                       ! If zeqn_2 < 0, then we adjust zh_min:
                                       ! we can be sure that zh_min < zh_2 < zh_max.
   zh_min = zh_2
                                       ! for zh_1, try 25% (0.1 pH units) above the current zh_min,
                                       ! but stay below SQRT(zh_min*zh_max) which would be equivalent
                                       ! to a bisection step on [pH@zh_min, pH@zh_max]
   zh_1   = MIN(zh_min*1.25_rk, SQRT(zh_min*zh_max))
   
ELSE ! we have got the root; unlikely, but one never knows

   SOLVE_AT_GENERAL_SEC = zh_2
   IF(PRESENT(p_val)) p_val = zeqn_2

   RETURN

ENDIF

! We now have the first pair completed (zh_2, zeqn_2).
! Define the second one (zh_1, zeqn_1), which is also the first iterate.
! zh_1 has already been set above

niter_atsec = 1                        ! Update counter of iterations


zeqn_1 = EQUATION_AT(p_alktot, zh_1,       p_dictot, p_bortot,                 &
                     p_po4tot, p_siltot, p_nh4tot, p_h2stot,                   &
                     p_so4tot, p_flutot)

! Adapt bracketing interval: we know that zh_1 <= zh <= zh_max (if zeqn_1 > 0)
! or zh_min <= zh <= zh_1 (if zeqn_1 < 0), so this can always be done

IF(zeqn_1 > 0._rk) THEN

   zh_min = zh_1

ELSEIF(zeqn_1 < 0._rk) THEN

   zh_max = zh_1

ELSE ! zh_1 is the root

   SOLVE_AT_GENERAL_SEC = zh_1
   IF(PRESENT(p_val)) p_val = zeqn_1

ENDIF


IF(ABS(zeqn_1) > zeqn_absmin) THEN     ! Swap zh_2 and zh_1 if ABS(zeqn_2) < ABS(zeqn_1)
                                       ! so that zh_2 and zh_1 lead to decreasing equation
                                       ! values (in absolute value)

   zh     = zh_1
   zeqn   = zeqn_1
   zh_1   = zh_2
   zeqn_1 = zeqn_2
   zh_2   = zh
   zeqn_2 = zeqn

ELSE

   zeqn_absmin = ABS(zeqn_1)

ENDIF


! Pre-calculate the first secant iterate (this is the second iterate)

niter_atsec = 2

zh_delta = -zeqn_1/((zeqn_2-zeqn_1)/(zh_2 - zh_1))

zh = zh_1 + zh_delta


! Make sure that zh_min < zh < zh_max (if not,
! bisect around zh_1 which is the best estimate)

IF (zh > zh_max) THEN                  ! this can only happen if zh_2 < zh_1
                                       ! and zeqn_2 > zeqn_1 > 0

   zh = SQRT(zh_1*zh_max)

ENDIF

IF (zh < zh_min) THEN                  ! this can only happen if zh_2 > zh_1
                                       ! and zeqn_2 < zeqn_1 < 0

   zh = SQRT(zh_1*zh_min)

ENDIF



DO

   IF(niter_atsec >= jp_maxniter_atsec) THEN
      zh = -1._rk
      EXIT
   ENDIF

   zeqn = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot, p_nh4tot, p_h2stot,                  &
                      p_so4tot, p_flutot)

   ! Adapt bracketing interval: since initially, zh_min <= zh <= zh_max
   ! we are sure that zh will improve either bracket, depending on the sign
   ! of zeqn
   IF(zeqn > 0._rk) THEN
     zh_min = zh
   ELSEIF(zeqn < 0._rk) THEN
     zh_max = zh
   ELSE
     ! zh is the root
     EXIT
   ENDIF


   ! start calculation of next iterate

   niter_atsec = niter_atsec + 1

   zh_2   = zh_1
   zeqn_2 = zeqn_1
   zh_1   = zh
   zeqn_1 = zeqn


   IF(ABS(zeqn) >= 0.5_rk*zeqn_absmin) THEN

      ! if the function evaluation at the current point
      ! is not decreasing faster in absolute value than
      ! we may expect for a bisection step, then take
      ! one bisection step on [ph_min, ph_max]
      ! ph_new = (ph_min + ph_max)/2d0
      ! In terms of [H]_new:
      ! [H]_new = 10**(-ph_new)
      !         = 10**(-(ph_min + ph_max)/2d0)
      !         = SQRT(10**(-(ph_min + phmax)))
      !         = SQRT(zh_max * zh_min)
     
      zh                = SQRT(zh_max * zh_min)
      zh_delta           = zh - zh_1

   ELSE


      ! \Delta H = -zeqn_1*(h_2 - h_1)/(zeqn_2 - zeqn_1) 
      ! H_new = H_1 + \Delta H

      zh_delta = -zeqn_1/((zeqn_2-zeqn_1)/(zh_2 - zh_1))

      zh       = zh_1 + zh_delta

#if defined(DEBUG_PHSOLVERS)
      PRINT*, '[SOLVE_AT_GENERAL_SEC] testing zh :', zh, zeqn, zh_delta
#endif


      IF( zh < zh_min ) THEN
         ! if [H]_new < [H]_min
         ! i.e., if ph_new > ph_max then
         ! take one bisection step on [ph_prev, ph_max]
         ! ph_new = (ph_prev + ph_max)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_max)/2d0)
         !         = SQRT(10**(-(ph_prev + phmax)))
         !         = SQRT([H]_old*10**(-ph_max))
         !         = SQRT([H]_old * zh_min)

         zh                = SQRT(zh_1 * zh_min)
         zh_delta          = zh - zh_1

      ENDIF

      IF( zh > zh_max ) THEN
         ! if [H]_new > [H]_max
         ! i.e., if ph_new < ph_min, then
         ! take one bisection step on [ph_min, ph_prev]
         ! ph_new = (ph_prev + ph_min)/2d0
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_prev + ph_min)/2d0)
         !         = SQRT(10**(-(ph_prev + ph_min)))
         !         = SQRT([H]_old*10**(-ph_min))
         !         = SQRT([H]_old * zhmax)

         zh                 = SQRT(zh_1 * zh_max)
         zh_delta           = zh - zh_1

      ENDIF

   ENDIF

   zeqn_absmin = MIN(ABS(zeqn), zeqn_absmin)


   ! Stop iterations once |([H]-[H_1])/[H_1]| < rdel

   l_exitnow = (ABS(zh_delta) < pp_rdel_ah_target*zh_1)

   IF(l_exitnow) EXIT

ENDDO


SOLVE_AT_GENERAL_SEC = zh

IF(PRESENT(p_val)) THEN

   IF(zh > 0._rk) THEN

      p_val = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot, p_nh4tot, p_h2stot,              &
                          p_so4tot, p_flutot)

   ELSE

     p_val = HUGE(1._rk)

   ENDIF

ENDIF


RETURN

!===============================================================================
 END FUNCTION SOLVE_AT_GENERAL_SEC
!===============================================================================





!===============================================================================
 FUNCTION SOLVE_AT_OCMIP(  p_alktot, p_dictot, p_bortot,                       &
                           p_po4tot, p_siltot, p_nh4tot, p_h2stot,             &
                           p_so4tot, p_flutot, p_hini,   p_val)
!===============================================================================

! Re-implementation and adaptation of the OCMIP-2 standard protocol solver.
! Requires a pair of bracketing values (interval) in the optional p_hini.

! Adaptations:
! * if p_hini is not provided, the extended cubic initialisation procedure
!   is used, and the interval set to that value +/- 0.5 pH units;
!
! * unlike the original OCMIP solver, an initial check regarding the validity
!   of the bracketing interval is carried out.
!
! * the convergence criterion was changed from the original absolute
!   error of 10^-10 (still present here in comments) to the same relative error
!   criterion used in all other solvers in this module (far more stringent).

! Returns
! * the value of the root ([H^+]) once the convergence criterion is verified
! * the value -1 if the bracketing interval is invalid or if convergence was
!   not detected (i.e., maximum number of iterations exceeded).

IMPLICIT NONE

REAL(KIND=rk) :: SOLVE_AT_OCMIP


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)            :: p_alktot
REAL(KIND=rk), INTENT(IN)            :: p_dictot
REAL(KIND=rk), INTENT(IN)            :: p_bortot
REAL(KIND=rk), INTENT(IN)            :: p_po4tot
REAL(KIND=rk), INTENT(IN)            :: p_siltot
REAL(KIND=rk), INTENT(IN)            :: p_nh4tot
REAL(KIND=rk), INTENT(IN)            :: p_h2stot
REAL(KIND=rk), INTENT(IN)            :: p_so4tot
REAL(KIND=rk), INTENT(IN)            :: p_flutot
REAL(KIND=rk), DIMENSION(2), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=rk), INTENT(OUT), OPTIONAL :: p_val


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk)  ::  zh_ini, zh, zh_prev, zh_lnfactor
REAL(KIND=rk)  ::  zh_low, zh_high
REAL(KIND=rk)  ::  zdelta, zh_delta, zh_delta_prev
REAL(KIND=rk)  ::  zeqn, zeqn_min, zeqn_max, zdeqndh

LOGICAL        ::  l_exitnow


!==============================================================================


IF(PRESENT(p_hini)) THEN

   zh_low  = MINVAL(p_hini(:))
   zh_high = MAXVAL(p_hini(:))
   zh_ini  = SQRT(zh_low*zh_high)      ! Average of pH@zh_low and pH@zh_high

ELSE

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_OCMIP] Calling AHINI_FOR_AT for h_ini'
#endif

   CALL AHINI_FOR_AT(p_alktot, p_dictot, p_bortot, zh_ini)

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_OCMIP] h_ini :', zh_ini
#endif

   zh_low  = zh_ini / 3.16_rk          ! +0.5 pH unit
   zh_high = zh_ini * 3.16_rk          ! -0.5 pH unit

ENDIF



niter_ocmip        = 0                 ! Reset counters of iterations


zeqn_max = EQUATION_AT(p_alktot, zh_low,   p_dictot, p_bortot,                 &
                       p_po4tot, p_siltot, p_nh4tot, p_h2stot,                 &
                       p_so4tot, p_flutot)


zeqn_min = EQUATION_AT(p_alktot, zh_high,  p_dictot, p_bortot,                 &
                       p_po4tot, p_siltot, p_nh4tot, p_h2stot,                 &
                       p_so4tot, p_flutot)

IF((zeqn_max < 0._rk) .OR. (zeqn_min > 0._rk)) THEN

   ! [zh_low, zh_high] does not bracket the root
   ! return to caller
   SOLVE_AT_OCMIP = -1._rk
   IF(PRESENT(p_val)) p_val = HUGE(1._rk)
   RETURN

ENDIF
      
zh            = 0.5_rk*(zh_low + zh_high)
zh_delta      = ABS(zh_high-zh_low)
zh_delta_prev = zh_delta

zeqn = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,                     &
                   p_po4tot, p_siltot, p_nh4tot, p_h2stot,                     &
                   p_so4tot, p_flutot, P_DERIVEQN = zdeqndh)

DO

   niter_ocmip = niter_ocmip + 1

   IF(niter_ocmip > jp_maxniter_ocmip) THEN ! too many iterations
      zh = -1._rk
      EXIT
   ENDIF

   zh_prev = zh

   IF(      ( ((zh-zh_low)*zdeqndh-zeqn)*((zh-zh_high)*zdeqndh-zeqn) > 0._rk ) &
       .OR. ( ABS(2.0*zeqn) > ABS(zh_delta_prev*zdeqndh)) ) THEN
      
      zh_delta_prev = zh_delta
      zh_delta      = 0.5_rk*(zh_high - zh_low)
      zh            = zh_low + zh_delta

      IF(zh_low == zh) EXIT ! zh_delta is not significant w/r to zh_low

   ELSE

      zh_delta_prev = zh_delta
      zh_delta      = zeqn/zdeqndh
      zh            = zh_prev - zh_delta

      IF(zh_prev == zh) EXIT ! zh_delta is not significant w/r to zh_prev

   END IF


   !l_exitnow = (ABS(zh_delta) < 1E-10_rk) ! original OCMIP has absolute criterion!
   l_exitnow = (ABS(zh_delta) < pp_rdel_ah_target*zh)

   IF(l_exitnow) EXIT

   zeqn = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot, p_nh4tot, p_h2stot,                  &
                      p_so4tot, p_flutot, P_DERIVEQN = zdeqndh)

   IF(zeqn < 0._rk) THEN

      zh_high  = zh

   ELSE

      zh_low   = zh

   END IF


ENDDO

SOLVE_AT_OCMIP = zh

IF(PRESENT(p_val)) THEN

   IF(zh > 0._rk) THEN

      p_val = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot, p_nh4tot, p_h2stot,              &
                          p_so4tot, p_flutot)

   ELSE

      p_val = HUGE(1._rk)

   ENDIF

ENDIF

RETURN

!===============================================================================
 END FUNCTION SOLVE_AT_OCMIP
!===============================================================================





!===============================================================================
 FUNCTION SOLVE_AT_FAST(p_alktot, p_dictot, p_bortot,                          &
                        p_po4tot, p_siltot, p_nh4tot, p_h2stot,                &
                        p_so4tot, p_flutot, p_hini,   p_val)
!===============================================================================

! Fast version of SOLVE_AT_GENERAL, without any bounds checking.

IMPLICIT NONE

REAL(KIND=rk) :: SOLVE_AT_FAST


!--------------------!
! Argument variables !
!--------------------!

REAL(KIND=rk), INTENT(IN)            :: p_alktot
REAL(KIND=rk), INTENT(IN)            :: p_dictot
REAL(KIND=rk), INTENT(IN)            :: p_bortot
REAL(KIND=rk), INTENT(IN)            :: p_po4tot
REAL(KIND=rk), INTENT(IN)            :: p_siltot
REAL(KIND=rk), INTENT(IN)            :: p_nh4tot
REAL(KIND=rk), INTENT(IN)            :: p_h2stot
REAL(KIND=rk), INTENT(IN)            :: p_so4tot
REAL(KIND=rk), INTENT(IN)            :: p_flutot
REAL(KIND=rk), INTENT(IN), OPTIONAL  :: p_hini
REAL(KIND=rk), INTENT(OUT), OPTIONAL :: p_val


!-----------------!
! Local variables !
!-----------------!

REAL(KIND=rk)  ::  zh_ini, zh, zh_prev, zh_lnfactor
REAL(KIND=rk)  ::  zhdelta
REAL(KIND=rk)  ::  zeqn, zdeqndh

LOGICAL        :: l_exitnow
REAL(KIND=rk), PARAMETER :: pz_exp_threshold = 1.0_rk


!==============================================================================


IF(PRESENT(p_hini)) THEN

   zh_ini = p_hini

ELSE

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_FAST] Calling AHINI_FOR_AT for h_ini'
#endif

   CALL AHINI_FOR_AT(p_alktot, p_dictot, p_bortot, zh_ini)

#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_FAST] h_ini :', zh_ini
#endif


ENDIF


zh = zh_ini

niter_atfast    = 0                 ! Reset counters of iterations


DO

   niter_atfast = niter_atfast + 1

   IF(niter_atfast > jp_maxniter_atfast) THEN
      zh = -1._rk
      EXIT
   ENDIF

   zh_prev = zh

   zeqn = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,                  &
                      p_po4tot, p_siltot, p_nh4tot, p_h2stot,                  &
                      p_so4tot, p_flutot, P_DERIVEQN = zdeqndh)

   IF(zeqn == 0._rk) EXIT               ! zh is the root

   zh_lnfactor = -zeqn/(zdeqndh*zh_prev)
   IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
      zh          = zh_prev*EXP(zh_lnfactor)
   ELSE
      zhdelta     = zh_lnfactor*zh_prev
      zh          = zh_prev + zhdelta
   ENDIF


#if defined(DEBUG_PHSOLVERS)
   PRINT*, '[SOLVE_AT_FAST] testing zh :', zh, zeqn, zh_lnfactor
#endif



   ! Stop iterations once |\delta{[H]}/[H]| < rdel
   ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
   ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

   ! Alternatively:
   ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
   !             ~ 1/LOG(10) * |\Delta [H]|/[H]
   !             < 1/LOG(10) * rdel

   ! Hence |zeqn/(zdeqndh*zh)| < rdel

   ! rdel <- pp_rdel_ah_target

   l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

   IF(l_exitnow) EXIT

ENDDO

SOLVE_AT_FAST = zh

IF(PRESENT(p_val)) THEN

   IF(zh > 0._rk) THEN

      p_val = EQUATION_AT(p_alktot, zh,       p_dictot, p_bortot,              &
                          p_po4tot, p_siltot, p_nh4tot, p_h2stot,              &
                          p_so4tot, p_flutot)

   ELSE

      p_val = HUGE(1._rk)

   ENDIF

ENDIF

RETURN

!===============================================================================
 END FUNCTION SOLVE_AT_FAST
!===============================================================================



END MODULE MOD_PHSOLVERS

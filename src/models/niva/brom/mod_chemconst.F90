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





MODULE MOD_CHEMCONST

use fabm_types, only:rk

IMPLICIT NONE

! --------------------------------------------------------
! List of subroutines for the chemical constants (PRIVATE)
! --------------------------------------------------------

PRIVATE AK_CARB_0_WEIS74
PRIVATE AK_CARB_1_MILL95, AK_CARB_2_MILL95
PRIVATE AK_CARB_1_LUEK00, AK_CARB_2_LUEK00
PRIVATE AK_CARB_1_ROYE93, AK_CARB_2_ROYE93
PRIVATE AK_BORA_DICK90
PRIVATE AK_PHOS_1_MILL95, AK_PHOS_2_MILL95, AK_PHOS_3_MILL95
PRIVATE AK_SILI_1_MILL95
PRIVATE AK_H2S_1_MILL95
PRIVATE AK_AMMO_1_YAMI95
PRIVATE AK_W_MILL95
PRIVATE AK_HSO4_DICK90
PRIVATE ABETA_HF_DIRI79
PRIVATE AK_HF_PEFR87


! --------------------------------------
! Parameters for usage within the module
! --------------------------------------

! Gas constant
! ------------

REAL(KIND=rk), PARAMETER, PRIVATE :: gasconst_bar_cm3_o_mol_k = 83.14510_rk ! libthdyct
!REAL(KIND=rk), PARAMETER, PRIVATE :: gasconst_bar_cm3_o_mol_k = 83.14472_rk ! Handbook (2007)

! 0 degrees centigrade in Kelvin
! ------------------------------

REAL(KIND=rk), PARAMETER, PRIVATE :: t_k_zerodegc = 273.15_rk ! Handbook (2007)


! --------------------------------------------------------------
! Chemical constants' products: for usage by users of the module
! --------------------------------------------------------------

! For each acid system A, 
! - api1_aaa <-- K_A1
! - api2_aaa <-- K_A1*K_A2
! - api3_aaa <-- K_A1*K_A2*K_A3
! - ...

REAL(KIND=rk) :: api1_dic, api2_dic
REAL(KIND=rk) :: api1_bor
REAL(KIND=rk) :: api1_po4, api2_po4, api3_po4
REAL(KIND=rk) :: api1_sil
REAL(KIND=rk) :: api1_nh4
REAL(KIND=rk) :: api1_h2s
REAL(KIND=rk) :: api1_so4
REAL(KIND=rk) :: api1_flu
REAL(KIND=rk) :: api1_wat, aphscale


!*******************************************************************************
CONTAINS
!*******************************************************************************

!=======================================================================
 SUBROUTINE SETUP_API4PHTOT(t_k, s, p_bar, &
            kc1, kc2, kb, kp1, kp2, kp3, ksi, & 
            knh4, kh2s, kw, kc0)
!=======================================================================

IMPLICIT NONE

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in Kelvin
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s
REAL(KIND=rk), INTENT(IN) :: p_bar

real(rk), optional, intent(out):: kc0, kc1, kc2, kb
real(rk), optional, intent(out):: kp1, kp2, kp3, ksi
real(rk), optional, intent(out):: knh4, kh2s, kw

! ---------------
! Local variables
! ---------------

REAL(KIND=rk) :: zcvt_htot_o_hsws, zcvt_htot_o_hfree


zcvt_htot_o_hsws  = 1._rk/ACVT_HSWS_O_HTOT(t_k, s, p_bar)
zcvt_htot_o_hfree = ACVT_HTOT_O_HFREE(t_k, s, p_bar)


api1_dic =            AK_CARB_1_LUEK00(t_k, s, p_bar)
kc1 = api1_dic

api2_dic = api1_dic * AK_CARB_2_LUEK00(t_k, s, p_bar)
kc2 = api2_dic

api1_bor =            AK_BORA_DICK90(t_k, s, p_bar)
kb = api1_bor

api1_po4 =            AK_PHOS_1_MILL95(t_k, s, p_bar) * zcvt_htot_o_hsws
kp1 = api1_po4
api2_po4 = api1_po4 * AK_PHOS_2_MILL95(t_k, s, p_bar) * zcvt_htot_o_hsws
kp2 = api2_po4
api3_po4 = api2_po4 * AK_PHOS_3_MILL95(t_k, s, p_bar) * zcvt_htot_o_hsws
kp3 = api3_po4

api1_sil =            AK_SILI_1_MILL95(t_k, s       ) * zcvt_htot_o_hsws
ksi = api1_sil

api1_nh4 =            AK_AMMO_1_YAMI95(t_k, s, p_bar) * zcvt_htot_o_hsws
knh4 = api1_nh4

api1_h2s =            AK_H2S_1_MILL95 (t_k, s, p_bar) * zcvt_htot_o_hsws
kh2s = api1_h2s

api1_so4 =            AK_HSO4_DICK90(t_k, s, p_bar) * zcvt_htot_o_hfree

api1_flu =            AK_HF_PEFR87(t_k, s, p_bar)

api1_wat =            AK_W_MILL95(t_k, s, p_bar)      * zcvt_htot_o_hsws
kw = api1_wat

aphscale =            zcvt_htot_o_hfree

kc0 = AK_CARB_0_WEIS74(t_k, s)

!=======================================================================
 END SUBROUTINE SETUP_API4PHTOT
!=======================================================================




!=======================================================================
 SUBROUTINE SETUP_API4PHSWS(t_k, s, p_bar)
!=======================================================================

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in Kelvin
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s
REAL(KIND=rk), INTENT(IN) :: p_bar

! ---------------
! Local variables
! ---------------

REAL(KIND=rk) :: zcvt_hsws_o_htot, zcvt_hsws_o_hfree


zcvt_hsws_o_htot  = ACVT_HSWS_O_HTOT(t_k, s, p_bar)
zcvt_hsws_o_hfree = ACVT_HSWS_O_HFREE(t_k, s, p_bar)


api1_dic =            AK_CARB_1_MILL95(t_k, s, p_bar)
api2_dic = api1_dic * AK_CARB_2_MILL95(t_k, s, p_bar)

api1_bor =            AK_BORA_DICK90(t_k, s, p_bar) * zcvt_hsws_o_htot

api1_po4 =            AK_PHOS_1_MILL95(t_k, s, p_bar)
api2_po4 = api1_po4 * AK_PHOS_2_MILL95(t_k, s, p_bar)
api3_po4 = api2_po4 * AK_PHOS_3_MILL95(t_k, s, p_bar)

api1_sil =            AK_SILI_1_MILL95(t_k, s       )

api1_nh4 =            AK_AMMO_1_YAMI95(t_k, s, p_bar)

api1_h2s =            AK_H2S_1_MILL95 (t_k, s, p_bar)

api1_so4 =            AK_HSO4_DICK90(t_k, s, p_bar) * zcvt_hsws_o_hfree

api1_flu =            zcvt_hsws_o_hfree/ABETA_HF_DIRI79(t_k, s, p_bar)

api1_wat =            AK_W_MILL95(t_k, s, p_bar)

aphscale =            zcvt_hsws_o_hfree

!=======================================================================
 END SUBROUTINE SETUP_API4PHSWS
!=======================================================================




!=======================================================================
 FUNCTION AK_CARB_0_WEIS74(t_k, s)
!=======================================================================

! Function calculates K0 in (mol/kg-SW)/atmosphere

! References: Weiss (1979) [(mol/kg-SW)/atm]
! pH scale  : N/A
! Note      : currently no pressure correction


IMPLICIT NONE

REAL(KIND=rk) :: AK_CARB_0_WEIS74


! ------------------
! Argument variables
! ------------------

!     s      : salinity
!     t_k    : temperature in K

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s


! ---------------
! Local variables
! ---------------

!     zt_k_o_100   : zt_k/100

REAL(KIND=rk) :: zt_k_o_100



zt_k_o_100 = t_k/100._rk

AK_CARB_0_WEIS74                                                       &
          = EXP( -60.2409_rk + 93.4517_rk/zt_k_o_100                   &
                + 23.3585_rk*LOG(zt_k_o_100)                           &
                + (   0.023517_rk - 0.023656_rk*zt_k_o_100             &
                   + 0.0047036_rk*zt_k_o_100*zt_k_o_100)*s )


RETURN

!=======================================================================
 END FUNCTION AK_CARB_0_WEIS74
!=======================================================================




!=======================================================================
 FUNCTION AK_CARB_1_MILL95(t_k, s, p_bar)
!=======================================================================
! Function calculates first dissociation constant of carbonic acid
! in mol/kg-SW on the SWS pH-scale.

! References: Millero (1995, eq 50 -- ln K1(COM))
!             Millero (1982) pressure correction
! pH scale:   SWS


IMPLICIT NONE

REAL(KIND=rk) :: AK_CARB_1_MILL95


! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in Kelvin
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s
REAL(KIND=rk), INTENT(IN) :: p_bar


! ---------------
! Local variables
! ---------------

!     zrt            : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc        : temperature in degrees Celsius
!     zdvi           : volume change for ionization
!     zdki           : compressibility change for ionization
!     zsqrts         : square root of salinity
!     zds            : salinity-34.8
!     zln_kc1_p0     : ln(K_C1) at p_bar = 0
!     zln_kc1_pp     : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
REAL(KIND=rk) :: zln_kc1_p0, zln_kc1_pp


! ln(K_C1) value at p_bar = 0

zsqrts     = SQRT(s)

zln_kc1_p0 =     2.18867_rk                                                    &
             - 2275.0360_rk/t_k                                                &
             - 1.468591_rk*LOG(t_k)                                            &
             + ( -0.138681_rk - 9.33291_rk/t_k)*zsqrts                         &
             + 0.0726483_rk*s - 0.00574938_rk*s*zsqrts


! Pressure correction

zt_degc    = t_k - t_k_zerodegc
zds        = s - 34.8_rk
zrt        = gasconst_bar_cm3_o_mol_k * t_k

zdvi       =  -25.50_rk - 0.151_rk*zds + 0.1271_rk*zt_degc
zdki       = ( -3.08_rk - 0.578_rk*zds + 0.0877_rk*zt_degc)*1.0E-03_rk

zln_kc1_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_C1 value

AK_CARB_1_MILL95 = EXP( zln_kc1_p0 + zln_kc1_pp )

RETURN

!=======================================================================
 END FUNCTION AK_CARB_1_MILL95
!=======================================================================




!=======================================================================
 FUNCTION AK_CARB_2_MILL95(t_k, s, p_bar)
!=======================================================================
! Function calculates second dissociation constant K1
! in mol/kg-SW on the SWS pH-scale.

! References: Millero (1995, eq 51 -- ln K2(COM))
!             Millero (1979) pressure correction
! pH scale:   SWS


IMPLICIT NONE

REAL(KIND=rk) :: AK_CARB_2_MILL95


! Argument variables
! ------------------

!     t_k    : temperature in Kelvin
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s
REAL(KIND=rk), INTENT(IN) :: p_bar


! Local variables
! ---------------

!     zrt            : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc        : temperature in degrees Celsius
!     zdvi           : volume change for ionization
!     zdki           : compressibility change for ionization
!     zsqrts         : square root of salinity
!     zds            : salinity-34.8
!     zln_kc2_p0     : ln(K_C2) at p_bar = 0
!     zln_kc2_pp     : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
REAL(KIND=rk) :: zln_kc2_p0, zln_kc2_pp


! ln(K_C2) value at p_bar = 0

zsqrts     = SQRT(s)

zln_kc2_p0 =    -0.84226_rk                                                    &
             - 3741.1288_rk/t_k                                                &
             -  1.437139_rk*LOG(t_k)                                           &
             + (-0.128417_rk - 24.41239_rk/t_k)*zsqrts                         &
             +  0.1195308_rk*s                                                 &
             - 0.00912840_rk*s*zsqrts


! Pressure correction

zt_degc    = t_k - t_k_zerodegc
zds        = s - 34.8_rk
zrt        = gasconst_bar_cm3_o_mol_k * t_k

zdvi       =  -15.82_rk + 0.321_rk*zds - 0.0219_rk*zt_degc
zdki       =  ( 1.13_rk - 0.314_rk*zds - 0.1475_rk*zt_degc)*1.0E-03_rk

zln_kc2_pp =  (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_C2 value

AK_CARB_2_MILL95  = EXP( zln_kc2_p0 + zln_kc2_pp )

RETURN

!=======================================================================
 END FUNCTION AK_CARB_2_MILL95
!=======================================================================




!=======================================================================
 FUNCTION AK_CARB_1_LUEK00(t_k, s, p_bar)
!=======================================================================
! Function calculates first dissociation constant of carbonic acid
! in mol/kg-SW on the Total pH-scale.

! References: Luecker et al. (2000) -- also Handbook (2007)
!             Millero (1979) pressure correction
! pH scale:   Total


IMPLICIT NONE

REAL(KIND=rk) :: AK_CARB_1_LUEK00


! Argument variables
! ------------------

!     t_k    : temperature in Kelvin
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s
REAL(KIND=rk), INTENT(IN) :: p_bar


! Local variables
! ---------------

!     zrt            : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc        : temperature in degrees Celsius
!     zdvi           : volume change for ionization
!     zdki           : compressibility change for ionization
!     zds            : salinity-34.8
!     zlog10_kc1_p0  : log_10(k_C1) at p_bar = 0
!     zln_kc1_pp     : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
REAL(KIND=rk) :: zlog10_kc1_p0, zln_kc1_pp


! log_10(K_C1) value at p_bar = 0

zlog10_kc1_p0 =   61.2172_rk                                                   &
                - 3633.86_rk/t_k                                               &
                - 9.67770_rk*LOG(t_k)                                          &
                + s*(0.011555 - s*0.0001152_rk)


! Pressure correction

zt_degc    = t_k - t_k_zerodegc
zds        = s - 34.8_rk
zrt        = gasconst_bar_cm3_o_mol_k * t_k

zdvi       =  -25.50_rk - 0.151_rk*zds + 0.1271_rk*zt_degc
zdki       = ( -3.08_rk - 0.578_rk*zds + 0.0877_rk*zt_degc)*1.0E-03_rk

zln_kc1_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_C1 value

AK_CARB_1_LUEK00 = 10._rk**zlog10_kc1_p0 * EXP(zln_kc1_pp)

RETURN

!=======================================================================
 END FUNCTION AK_CARB_1_LUEK00
!=======================================================================




!=======================================================================
 FUNCTION AK_CARB_2_LUEK00(t_k, s, p_bar)
!=======================================================================
! Function calculates second dissociation constant K1
! in mol/kg-SW on the Total pH-scale.

! References: Luecker et al. (2000) -- also Handbook (2007)
!             Millero (1979) pressure correction
! pH scale:   Total


IMPLICIT NONE

REAL(KIND=rk) :: AK_CARB_2_LUEK00


! Argument variables
! ------------------

!     t_k    : temperature in Kelvin
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s
REAL(KIND=rk), INTENT(IN) :: p_bar


! Local variables
! ---------------

!     zrt            : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc        : temperature in degrees Celsius
!     zdvi           : volume change for ionization
!     zdki           : compressibility change for ionization
!     zsqrts         : square root of salinity
!     zds            : salinity-34.8
!     zlog10_kc2_p0  : log_10(K_C2) at p_bar = 0
!     zln_kc2_pp     : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
REAL(KIND=rk) :: zlog10_kc2_p0, zln_kc2_pp


! log_10(K_C2) value at p_bar = 0

zlog10_kc2_p0 =  -25.9290_rk                                                   &
                -  471.78_rk/t_k + 3.16967_rk*LOG(t_k)                         &
                + s*(0.01781_rk - s*0.0001122_rk)


! Pressure correction

zt_degc    = t_k - t_k_zerodegc
zds        = s - 34.8_rk
zrt        = gasconst_bar_cm3_o_mol_k * t_k

zdvi       =  -15.82_rk + 0.321_rk*zds - 0.0219_rk*zt_degc
zdki       =  ( 1.13_rk - 0.314_rk*zds - 0.1475_rk*zt_degc)*1.0E-03_rk

zln_kc2_pp =  (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_C2 value

AK_CARB_2_LUEK00  = 10._rk**zlog10_kc2_p0 *EXP(zln_kc2_pp)

RETURN

!=======================================================================
 END FUNCTION AK_CARB_2_LUEK00
!=======================================================================








!=======================================================================
 FUNCTION AK_CARB_1_ROYE93(t_k, s, p_bar)
!=======================================================================
! Function calculates first dissociation constant of carbonic acid
! in mol/kg-SW on the Total pH-scale.

! References: Roy et al. (1993) -- also Handbook (1994)
!             Millero (1979) pressure correction
! pH scale  : Total
! Note      : converted here from mol/kg-H2O to mol/kg-SW


IMPLICIT NONE

REAL(KIND=rk) :: AK_CARB_1_ROYE93


! Argument variables
! ------------------

!     t_k    : temperature in Kelvin
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s
REAL(KIND=rk), INTENT(IN) :: p_bar


! Local variables
! ---------------

!     zrt            : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc        : temperature in degrees Celsius
!     zdvi           : volume change for ionization
!     zdki           : compressibility change for ionization
!     zds            : salinity-34.8
!     zln_kc1_p0     : ln(k_C1) at p_bar = 0
!     zln_kc1_pp     : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zsqrts, zcvt_to_kgsw
REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds
REAL(KIND=rk) :: zln_kc1_p0, zln_kc1_pp


! ln(K_C1) value at p_bar = 0

zsqrts     = SQRT(s)
zcvt_to_kgsw    = ACVT_KGH2O_O_KGSW(s)

zln_kc1_p0 =   -2307.1255_rk/t_k + 2.83655_rk - 1.5529413_rk*LOG(t_k)        &
             + (-4.0484_rk/t_k - 0.20760841)*zsqrts                          &
             + 0.08468345*s                                                  &
             - 0.00654208*zsqrts*s


! Pressure correction

zt_degc    = t_k - t_k_zerodegc
zds        = s - 34.8_rk
zrt        = gasconst_bar_cm3_o_mol_k * t_k

zdvi       =  -25.50_rk - 0.151_rk*zds + 0.1271_rk*zt_degc
zdki       = ( -3.08_rk - 0.578_rk*zds + 0.0877_rk*zt_degc)*1.0E-03_rk

zln_kc1_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_C1 value

AK_CARB_1_ROYE93 = EXP(zln_kc1_p0 + zln_kc1_pp) * zcvt_to_kgsw

RETURN

!=======================================================================
 END FUNCTION AK_CARB_1_ROYE93
!=======================================================================




!=======================================================================
 FUNCTION AK_CARB_2_ROYE93(t_k, s, p_bar)
!=======================================================================
! Function calculates second dissociation constant K1
! in mol/kg-SW on the Total pH-scale.

! References: Roy et al. (1993) -- also Handbook (1994)
!             Millero (1979) pressure correction
! pH scale  : Total
! Note      : converted here from mol/kg-H2O to mol/kg-SW


IMPLICIT NONE

REAL(KIND=rk) :: AK_CARB_2_ROYE93


! Argument variables
! ------------------

!     t_k    : temperature in Kelvin
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s
REAL(KIND=rk), INTENT(IN) :: p_bar


! Local variables
! ---------------

!     zrt            : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc        : temperature in degrees Celsius
!     zdvi           : volume change for ionization
!     zdki           : compressibility change for ionization
!     zsqrts         : square root of salinity
!     zds            : salinity-34.8
!     zln_kc2_p0     : ln(K_C2) at p_bar = 0
!     zln_kc2_pp     : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zsqrts, zcvt_to_kgsw
REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds
REAL(KIND=rk) :: zln_kc2_p0, zln_kc2_pp


! ln(K_C2) value at p_bar = 0

zsqrts     = SQRT(s)
zcvt_to_kgsw    = ACVT_KGH2O_O_KGSW(s)

zln_kc2_p0 =   -3351.6106_rk/t_k - 9.226508_rk - 0.2005743_rk*LOG(t_k)         &
             + ( -23.9722_rk/t_k - 0.106901773_rk)*zsqrts                      &
             + 0.1130822*s - 0.00846934_rk*zsqrts*s


! Pressure correction

zt_degc    = t_k - t_k_zerodegc
zds        = s - 34.8_rk
zrt        = gasconst_bar_cm3_o_mol_k * t_k

zdvi       =  -15.82_rk + 0.321_rk*zds - 0.0219_rk*zt_degc
zdki       =  ( 1.13_rk - 0.314_rk*zds - 0.1475_rk*zt_degc)*1.0E-03_rk

zln_kc2_pp =  (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_C2 value

AK_CARB_2_ROYE93  = EXP(zln_kc2_p0 + zln_kc2_pp) * zcvt_to_kgsw

RETURN

!=======================================================================
 END FUNCTION AK_CARB_2_ROYE93
!=======================================================================




!=======================================================================
 FUNCTION AK_BORA_DICK90(t_k, s, p_bar)
!=======================================================================
! Function calculates boric acid dissociation constant KB
! in mol/kg-SW on the total pH-scale.

! References: Dickson (1990, eq. 23) -- also Handbook (2007, eq. 37)
!             Millero (1979) pressure correction
! pH scale  : total


IMPLICIT NONE

REAL(KIND=rk) :: AK_BORA_DICK90


! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in Kelvin
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) :: t_k
REAL(KIND=rk), INTENT(IN) :: s
REAL(KIND=rk), INTENT(IN) :: p_bar


! ---------------
! Local variables
! ---------------

!     zrt          : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc      : temperature in degrees Celsius
!     zdvi         : volume change for ionization
!     zdki         : compressibility change for ionization
!     zsqrts       : square root of salinity
!     zds          : salinity-34.8
!     zln_kb_p0    : K_b at p_bar = 0
!     zln_kb_pp    : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
REAL(KIND=rk) :: zln_kb_p0, zln_kb_pp


! ln(K_B) value at p_bar = 0

zsqrts     = SQRT(s)

zln_kb_p0  = ( -8966.90_rk                                                     &
                          + zsqrts*( -2890.53_rk                               &
                          + zsqrts*(  -77.942_rk                               &
                          + zsqrts*(    1.728_rk - 0.0996_rk*zsqrts)))) / t_k  &
               +  148.0248_rk + zsqrts*(137.1942_rk + zsqrts*1.62142_rk)       &
               + (-24.4344_rk + zsqrts*(-25.085_rk - zsqrts*0.2474_rk)) * LOG(t_k) &
               + 0.053105_rk*zsqrts*t_k


! Pressure correction

zt_degc   = t_k - t_k_zerodegc
zds       = s - 34.8_rk
zrt       = gasconst_bar_cm3_o_mol_k * t_k

zdvi      = -29.48_rk + 0.295_rk*zds + 0.1622_rk*zt_degc - 0.002608_rk*zt_degc*zt_degc
zdki      = (-2.84_rk + 0.354_rk*zds)*1.0E-03_rk

zln_kb_pp =  (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_B value

AK_BORA_DICK90   = EXP( zln_kb_p0 + zln_kb_pp )

!=======================================================================
 END FUNCTION AK_BORA_DICK90
!=======================================================================




!=======================================================================
 FUNCTION AK_W_MILL95(t_k, s, p_bar)
!=======================================================================

! Function calculates water dissociation constant Kw in (mol/kg-SW)^2

! References: Millero (1995) for value at p_bar = 0
!             Millero (pers. comm. 1996) for pressure correction
! pH scale  : SWS


IMPLICIT NONE

REAL(KIND=rk) :: AK_W_MILL95


! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     zrt        : R*t_k
!     zt_degc    : temperature in degrees Celsius
!     zdvi       : volume change for ionization
!     zdki       : compressibility change for ionization
!     zln_kw_p0  : ln(K_w) at p_bar = 0
!     zln_kw_pp  : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
REAL(KIND=rk) :: zln_kw_p0, zln_kw_pp


! ln(K_w) value at p_bar = 0

zln_kw_p0 =    148.9802_rk                                                     &
             - 13847.26_rk/t_k                                                 &
             -  23.6521_rk*LOG(t_k)                                            &
             + ( -5.977_rk + 118.67_rk/t_k + 1.0495_rk*LOG(t_k))*SQRT(s)       &
             - 0.01615_rk*s


! Pressure correction

zt_degc = t_k - t_k_zerodegc
zrt     = gasconst_bar_cm3_o_mol_k * t_k

zdvi    =  -20.02_rk + 0.1119_rk*zt_degc - 0.1409E-02_rk*zt_degc*zt_degc
zdki    = ( -5.13_rk + 0.0794_rk*zt_degc)*1.0E-03_rk

zln_kw_pp =  (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_w value

AK_W_MILL95 = EXP( zln_kw_p0 + zln_kw_pp )


RETURN

!=======================================================================
END FUNCTION AK_W_MILL95
!=======================================================================




!=======================================================================
FUNCTION AK_PHOS_1_MILL95(t_k, s, p_bar)
!=======================================================================

! Function returns the first dissociation constant
! of phosphoric acid (H3PO4) in seawater

! References: Yao and Millero (1995)
!             Millero (1995) for pressure correction
! pH scale  : SWS


IMPLICIT NONE

REAL(KIND=rk) :: AK_PHOS_1_MILL95

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar

! ---------------
! Local variables
! ---------------

!     zrt          : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc      : temperature in degrees Celsius
!     zdvi         : volume change for ionization
!     zdki         : compressibility change for ionization
!     zln_kp1_p0   : ln(K_p1) at p_bar = 0
!     zln_kp1_pp   : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki
REAL(KIND=rk) :: zln_kp1_p0, zln_kp1_pp


! ln(K_P1) for p_bar = 0

zln_kp1_p0 =      115.54_rk - 4576.752_rk/t_k - 18.453_rk*LOG(t_k)     &
             + ( 0.69171_rk -  106.736_rk/t_k)* SQRT(s)               &
             + (-0.01844_rk -  0.65643_rk/t_k)*s 


! Pressure correction

zt_degc   = t_k - t_k_zerodegc
zrt       = gasconst_bar_cm3_o_mol_k * t_k

zdvi      =  -14.51_rk + 0.1211_rk*zt_degc - 0.321E-03*zt_degc*zt_degc
zdki      = ( -2.67_rk + 0.0427_rk*zt_degc)*1.0E-03_rk

zln_kp1_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final value of K_P1

AK_PHOS_1_MILL95 = EXP(zln_kp1_p0 + zln_kp1_pp)

RETURN

!=======================================================================
 END FUNCTION AK_PHOS_1_MILL95
!=======================================================================



!=======================================================================
 FUNCTION AK_PHOS_2_MILL95(t_k, s, p_bar)
!=======================================================================

! Function returns the second dissociation constant
! of phosphoric acid (H3PO4) in seawater

! References: Yao and Millero (1995)
!             Millero (1995) for pressure correction
! pH scale  : SWS


IMPLICIT NONE

REAL(KIND=rk) :: AK_PHOS_2_MILL95

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     zrt          : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc      : temperature in degrees Celsius
!     zdvi         : volume change for ionization
!     zdki         : compressibility change for ionization
!     zln_kp2_p0   : ln(K_P2) at p_bar = 0
!     zln_kp2_pp   : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki
REAL(KIND=rk) :: zln_kp2_p0, zln_kp2_pp


! ln(K_P2) for p_bar = 0

zln_kp2_p0 =   172.1033_rk                                                     &
             - 8814.715_rk/t_k                                                 &
             -   27.927_rk*LOG(t_k)                                            &
             + (  1.3566_rk -  160.340_rk/t_k)*SQRT(s)                         &
             + (-0.05778_rk +  0.37335_rk/t_k)*s


! Pressure correction

zt_degc    = t_k - t_k_zerodegc
zrt        = gasconst_bar_cm3_o_mol_k * t_k

zdvi       =  -23.12_rk + 0.1758_rk*zt_degc -2.647E-03_rk*zt_degc*zt_degc
zdki       = ( -5.15_rk +   0.09_rk*zt_degc)*1.0E-03_rk

zln_kp2_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_P2 value

AK_PHOS_2_MILL95  = EXP( zln_kp2_p0 + zln_kp2_pp )

RETURN

!=======================================================================
 END FUNCTION AK_PHOS_2_MILL95
!=======================================================================



!=======================================================================
FUNCTION AK_PHOS_3_MILL95(t_k, s, p_bar)
!=======================================================================

! Function returns the third dissociation constant
! of phosphoric acid (H3PO4) in seawater

! References: Yao and Millero (1995)
!             Millero (1995) for pressure correction
! pH scale  : SWS


IMPLICIT NONE

REAL(KIND=rk) :: AK_PHOS_3_MILL95

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     zrt          : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc      : temperature in degrees Celsius
!     zdvi         : volume change for ionization
!     zdki         : compressibility change for ionization
!     zln_kp3_p0   : ln(K_P3) at p_bar = 0
!     zln_kp3_pp   : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki
REAL(KIND=rk) :: zln_kp3_p0, zln_kp3_pp


! ln(K_P3) for p_bar = 0

zln_kp3_p0 =     -18.126_rk  -  3070.75_rk/t_k                        &
               + ( 2.81197_rk + 17.27039_rk/t_k)*SQRT(s)               &
               + (-0.09984_rk - 44.99486_rk/t_k)*s


! Pressure correction

zt_degc   = t_k - t_k_zerodegc
zrt       = gasconst_bar_cm3_o_mol_k * t_k

zdvi      =  -26.57_rk + 0.2020_rk*zt_degc -3.042E-03*zt_degc*zt_degc
zdki      = ( -4.08_rk + 0.0714_rk*zt_degc)*1.0E-03_rk

zln_kp3_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_P3 value

AK_PHOS_3_MILL95 = EXP( zln_kp3_p0 + zln_kp3_pp )

RETURN

!=======================================================================
END FUNCTION AK_PHOS_3_MILL95
!=======================================================================



!=======================================================================
FUNCTION AK_SILI_1_MILL95(t_k, s)
!=======================================================================

! Function returns the first dissociation constant
! of silicic acid (H4SiO4) in seawater

! References: Yao and Millero (1995) cited by Millero (1995)
! pH scale  : SWS (according to Dickson et al, 2007)
! Note      : No pressure correction available
! Note      : converted here from mol/kg-H2O to mol/kg-sw


IMPLICIT NONE

REAL(KIND=rk) :: AK_SILI_1_MILL95

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity

REAL(KIND=rk), INTENT(IN) :: t_k, s


! ---------------
! Local variables
! ---------------

!     zcvt_to_kgsw: fraction of pure water in 1 kg seawater at salinity s
!     zionst      : ionic strength [mol/kg-H2O]
!     zln_ksi1_p0 : ln(K_Si1) at p_bar = 0
!     zln_ksi1_pp : pressure correciotn for p_bar /= 0

REAL(KIND=rk) :: zionst, zcvt_to_kgsw
REAL(KIND=rk) :: zln_ksi1_p0, zln_ksi1_pp


! K_Si1 value at p_bar = 0

zcvt_to_kgsw = ACVT_KGH2O_O_KGSW(s)
zionst       = A_IONSTRENGTH_SALIN(s)/zcvt_to_kgsw  ! mol/kg-H2O !!

zln_ksi1_p0  =     117.40_rk -  8904.2_rk/t_k - 19.334_rk * LOG(t_k)    &
               + ( 3.5913_rk -  458.79_rk/t_k) * SQRT(zionst)           &
               + (-1.5998_rk +  188.74_rk/t_k) * zionst                 &
               + (0.07871_rk - 12.1652_rk/t_k) * zionst*zionst 


! Pressure correction : currently none

zln_ksi1_pp = 0._rk


! Final value

AK_SILI_1_MILL95 = EXP( zln_ksi1_p0 + zln_ksi1_pp ) * zcvt_to_kgsw


RETURN

!=======================================================================
END FUNCTION AK_SILI_1_MILL95
!=======================================================================



!=======================================================================
FUNCTION AK_H2S_1_MILL95(t_k, s, p_bar)
!=======================================================================

! Function returns the dissociation constant of hydrogen sulfide in sea-water


! References: Millero et al. (1988) (cited by Millero (1995)
!             Millero (1995) for pressure correction
! pH scale  : - SWS (according to Yao and Millero, 1995, p. 82: "refitted if necessary")
!             - Total (according to Lewis and Wallace, 1998)
! Note      : we stick to SWS here for the time being
! Note      : the fits from Millero (1995) and Yao and Millero (1995)
!             derive from Millero et al. (1998), with all the coefficients
!             multiplied by -ln(10)

IMPLICIT NONE

REAL(KIND=rk) :: AK_H2S_1_MILL95

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     zt_degc      : temperature in degrees Celsius
!     zrt          : R*t_k, R in bar*cm3/(mol*K)
!     zdvi         : volume change for ionization
!     zdki         : compressibility change for ionization
!     zln_kh2s_p0  : ln(K_H2S) at p_bar = 0
!     zln_kh2s_pp  : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki
REAL(KIND=rk) :: zln_kh2s_p0, zln_kh2s_pp


! K_H2S value at p_bar = 0
! ------------------------

zln_kh2s_p0  =    225.838_rk                                                     &
                - 13275.3_rk/t_k                                                 &
                - 34.6435_rk * LOG(t_k)                                          &
                +  0.3449_rk*SQRT(s)                                             &
                -  0.0274_rk*s


! Pressure correction
! -------------------

zt_degc      = t_k - t_k_zerodegc
zrt          = gasconst_bar_cm3_o_mol_k * t_k

zdvi         =  -14.80_rk + zt_degc*(0.0020_rk - zt_degc*0.400E-03_rk)
zdki         = (  2.89_rk + zt_degc*0.054_rk)*1.0E-03_rk

zln_kh2s_pp  = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_H2S value
! -----------------

AK_H2S_1_MILL95 = EXP( zln_kh2s_p0 + zln_kh2s_pp )


RETURN

!=======================================================================
END FUNCTION AK_H2S_1_MILL95
!=======================================================================



!=======================================================================
FUNCTION AK_AMMO_1_YAMI95(t_k, s, p_bar)
!=======================================================================

! Function returns the dissociation constant
! of ammonium in sea-water [mol/kg-SW]


! References: Yao and Millero (1995)
!             Millero (1995) for pressure correction
! pH scale  : SWS


IMPLICIT NONE

REAL(KIND=rk) :: AK_AMMO_1_YAMI95

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     zt_degc      : temperature in degrees Celsius
!     zrt          : R*t_k, R in bar*cm3/(mol*K)
!     zdvi         : volume change for ionization
!     zdki         : compressibility change for ionization
!     zln_knh4_p0  : ln(K_NH4) at p_bar = 0
!     zln_knh4_pp  : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki
REAL(KIND=rk) :: zln_knh4_p0, zln_knh4_pp


! K_NH4 value at p_bar = 0
! ------------------------

zln_knh4_p0  =    -0.25444_rk -  6285.33_rk/t_k + 0.0001635_rk*t_k             &
               + ( 0.46532_rk - 123.7184_rk/t_k) * SQRT(s)                     &
               + (-0.01992_rk +  3.17556_rk/t_k) * s   


! Pressure correction
! -------------------

zt_degc      = t_k - t_k_zerodegc
zrt          = gasconst_bar_cm3_o_mol_k * t_k

zdvi         =  -26.43_rk + zt_degc*(0.0889_rk - zt_degc*0.905E-03_rk)
zdki         = ( -5.03_rk + zt_degc*0.0814_rk)*1.0E-03_rk

zln_knh4_pp  = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final K_NH4 value
! -----------------

AK_AMMO_1_YAMI95  = EXP( zln_knh4_p0 + zln_knh4_pp )

RETURN

!=======================================================================
END FUNCTION AK_AMMO_1_YAMI95
!=======================================================================



!=======================================================================
FUNCTION ACVT_KGH2O_O_KGSW(s)
!=======================================================================

! Function returns the mass of pure water in one kg of seawater
! of salinity s

! References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
!             "Handbook (2007)" -- Handbook (2007)
! pH scale:   N/A


IMPLICIT NONE

REAL(KIND=rk) :: ACVT_KGH2O_O_KGSW

REAL(KIND=rk), INTENT(IN) :: s

!ACVT_KGH2O_O_KGSW = 1._rk - 0.0010049_rk*s ! libthdyct
ACVT_KGH2O_O_KGSW = 1._rk - 0.001005_rk*s ! Handbook (2007)

RETURN


!=======================================================================
END FUNCTION ACVT_KGH2O_O_KGSW
!=======================================================================



!=======================================================================
FUNCTION A_IONSTRENGTH_SALIN(s)
!=======================================================================

! Function calculates ionic strength in mol/kg-SW, for given salinity.

! References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
!             "Handbook (2007)" -- Handbook (2007)
! pH scale:   N/A


IMPLICIT NONE

REAL(KIND=rk) :: A_IONSTRENGTH_SALIN


! ------------------
! Argument variables
! ------------------

REAL(KIND=rk), INTENT(IN) :: s



!A_IONSTRENGTH_SALIN = (0.019920D+00*s) ! libthdyct
A_IONSTRENGTH_SALIN = (0.019924D+00*s) ! Handbook (2007)

RETURN

!=======================================================================
END FUNCTION A_IONSTRENGTH_SALIN
!=======================================================================




!=======================================================================
FUNCTION ABETA_HF_DIRI79(t_k, s, p_bar)
!=======================================================================

! Function calculates association constant \beta_{HF} [(mol/kg-SW)^{-1}]
! in (mol/kg-SW)^{-1}, where
!   \beta_{HF} = \frac{ [HF] }{ [H^{+}] [F^{-}] }

! References: Dickson and Riley (1979)
!             Millero (1995) for pressure correction
! pH scale  : free
! Note      : converted here from mol/kg-H2O to mol/kg-SW

IMPLICIT NONE

REAL(KIND=rk) :: ABETA_HF_DIRI79


! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     zrt            : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc        : temperature in degrees Celsius
!     zdvi           : volume change for ionization
!     zdki           : compressibility change for ionization
!     zionst         : ionic strength [mol/kg-H2O]
!     zcvt_to_kgsw   : mass of pure water in 1kg of seawater as a fct. of salinity
!     zln_bhf_p0     : \beta_HF at p_bar = 0
!     zln_khf_pp     : pressure correction for k_HF = 1/\beta_HF at p_bar /= 0


REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
REAL(KIND=rk) :: zionst, zcvt_to_kgsw
REAL(KIND=rk) :: zln_bhf_p0, zln_khf_pp


! \beta_HF at p_bar = 0
! ---------------------

zcvt_to_kgsw    = ACVT_KGH2O_O_KGSW(s)
zionst          = A_IONSTRENGTH_SALIN(s)/zcvt_to_kgsw 

zln_bhf_p0      = -1590.2_rk/t_k + 12.641_rk - 1.525_rk*SQRT(zionst)


! Pressure correction
! -------------------

zt_degc      = t_k - t_k_zerodegc
zrt          = gasconst_bar_cm3_o_mol_k * t_k

zdvi         =   -9.78_rk + zt_degc*(-0.0090_rk - zt_degc*0.942E-03_rk)
zdki         = ( -3.91_rk + zt_degc*0.054_rk)*1.0E-03_rk

zln_khf_pp   = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final \beta_HF value
! --------------------
!  notice that  ln(k_HF(P)) = ln(k_HF(0)) + zln_khf_pp
!         <=>  -ln(\beta_HF(P)) = -ln(\beta_HF(0)) + zln_khf_pp
!         <=>   ln(\beta_HF(P)) =  ln(\beta_HF(0)) - zln_khf_pp

ABETA_HF_DIRI79 = EXP(zln_bhf_p0 - zln_khf_pp ) / zcvt_to_kgsw

RETURN

!=======================================================================
END FUNCTION ABETA_HF_DIRI79
!=======================================================================




!=======================================================================
FUNCTION AK_HF_PEFR87(t_k, s, p_bar)
!=======================================================================

! Function calculates dissociation constant for hydrogen fluoride
! in mol/kg-SW

! References: Perez and Fraga (1987)
!             Millero (1995) for pressure correction
! pH scale  : Total (according to Handbook, 2007)

IMPLICIT NONE

REAL(KIND=rk) :: AK_HF_PEFR87


! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     zrt            : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc        : temperature in degrees Celsius
!     zdvi           : volume change for ionization
!     zdki           : compressibility change for ionization
!     zln_khf_p0     : ln(K_HF) at p_bar = 0
!     zln_khf_pp     : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki, zds, zsqrts
REAL(KIND=rk) :: zln_khf_p0, zln_khf_pp


! ln(K_HF) at p_bar = 0

zln_khf_p0   = 874._rk/t_k - 9.68_rk + 0.111_rk*SQRT(s)


! Pressure correction

zt_degc      = t_k - t_k_zerodegc
zrt          = gasconst_bar_cm3_o_mol_k * t_k

zdvi         =   -9.78_rk + zt_degc*(-0.0090_rk - zt_degc*0.942E-03_rk)
zdki         = ( -3.91_rk + zt_degc*0.054_rk)*1.0E-03_rk

zln_khf_pp   = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final value of K_HF

AK_HF_PEFR87 = EXP( zln_khf_p0 + zln_khf_pp )

RETURN

!=======================================================================
END FUNCTION AK_HF_PEFR87
!=======================================================================




!=======================================================================
FUNCTION AK_HSO4_DICK90(t_k, s, p_bar)
!=======================================================================

! Function returns the dissociation constant of hydrogen sulfate (bisulfate)

! References: Dickson (1990) -- also Handbook (2007)
!             Millero (1995) for pressure correction
! pH scale  : free
! Note      : converted here from mol/kg-H2O to mol/kg-SW

IMPLICIT NONE

REAL(KIND=rk) :: AK_HSO4_DICK90


! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     zrt            : R*t_k, R in bar*cm3/(mol*K)
!     zt_degc        : temperature in degrees Celsius
!     zdvi           : volume change for ionization
!     zdki           : compressibility change for ionization
!     zionst         : ionic strength in mol/-kg-H2O
!     zsqrti         : square root og ion strength
!     zcvt_to_kgsw   : mass of pure water in 1kg of seawater as a fct. of salinity
!     zln_khso4_p0   : K_HSO4 at p_bar = 0
!     zln_khso4_pp   : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zt_degc, zdvi, zdki
REAL(KIND=rk) :: zcvt_to_kgsw, zionst, zsqrti
REAL(KIND=rk) :: zln_khso4_p0, zln_khso4_pp


! ln(K_HSO4) at p_bar = 0

zcvt_to_kgsw = ACVT_KGH2O_O_KGSW(s)
zionst       = A_IONSTRENGTH_SALIN(s)/zcvt_to_kgsw
zsqrti       = SQRT(zionst)

zln_khso4_p0 =     -4276.1_rk/t_k + 141.328_rk -  23.093_rk*LOG(t_k)           &
               + (-13856._rk/t_k  +  324.57_rk -  47.986_rk*LOG(t_k)) * zsqrti &
               + ( 35474._rk/t_k  -  771.54_rk + 114.723_rk*LOG(t_k)) * zionst &
               - (  2698._rk/t_k)*zsqrti * zionst                              &
               + (  1776._rk/t_k)*zionst*zionst


! Pressure correction

zt_degc      = t_k - t_k_zerodegc
zrt          = gasconst_bar_cm3_o_mol_k * t_k

zdvi         =  -18.03_rk + zt_degc*(0.0466_rk + zt_degc*0.316E-03_rk)
zdki         = ( -4.53_rk + zt_degc*0.0900_rk)*1.0E-03_rk

zln_khso4_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! ln(K_HSO4) at p_bar = 0

AK_HSO4_DICK90 = zcvt_to_kgsw * EXP( zln_khso4_p0 + zln_khso4_pp )

RETURN

!=======================================================================
END FUNCTION AK_HSO4_DICK90
!=======================================================================




!=======================================================================
FUNCTION ASP_CALC_MUCC83(t_k, s, p_bar)
!=======================================================================

! Function returns stoechiometric solubility product
! of calcite in seawater

! References: Mucci (1983)
!             Millero (1995) for pressure correction
! pH scale  : N/A
! Units     : (mol/kg-SW)^2


IMPLICIT NONE

REAL(KIND=rk) :: ASP_CALC_MUCC83

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar

! ---------------
! Local variables
! ---------------

!     zrt          : R*t_k, R in bar*cm3/(mol*K)
!     zsqrts       : square root of salinity
!     zt_degc      : temperature in degrees Celsius
!     zdvi         : volume change for ionization
!     zdki         : compressibility change for ionization
!     zln_kp1_p0   : ln(K_p1) at p_bar = 0
!     zln_kp1_pp   : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zsqrts, zt_degc, zdvi, zdki
REAL(KIND=rk) :: zlog10_kspcalc_p0, zln_kspcalc_pp


zsqrts    = SQRT(s)

! log10(Ksp_Calc) for p_bar = 0
zlog10_kspcalc_p0 = &
             -171.9065_rk - 0.077993_rk*t_k                           &
            + 2839.319_rk/t_k + 71.595_rk*LOG10(t_k)                  &
            + ( -0.77712_rk + 0.0028426*t_k + 178.34_rk/t_k)*zsqrts   &
            - 0.07711_rk*s + 0.0041249_rk*s*zsqrts


! Pressure correction
zt_degc   = t_k - t_k_zerodegc
zrt       = gasconst_bar_cm3_o_mol_k * t_k

zdvi      =  -48.76_rk + 0.5304_rk*zt_degc
zdki      = (-11.76_rk + 0.3692_rk*zt_degc)*1.0E-03_rk

zln_kspcalc_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final value of Ksp_Calc

ASP_CALC_MUCC83 = 10._rk**(zlog10_kspcalc_p0) * EXP(zln_kspcalc_pp)

RETURN

!=======================================================================
 END FUNCTION ASP_CALC_MUCC83
!=======================================================================




!=======================================================================
FUNCTION ASP_ARAG_MUCC83(t_k, s, p_bar)
!=======================================================================

! Function returns stoechiometric solubility product
! of aragonite in seawater

! References: Mucci (1983)
!             Millero (1979) for pressure correction
! pH scale  : N/A
! Units     : (mol/kg-SW)^2


IMPLICIT NONE

REAL(KIND=rk) :: ASP_ARAG_MUCC83

! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar

! ---------------
! Local variables
! ---------------

!     zrt          : R*t_k, R in bar*cm3/(mol*K)
!     zsqrts       : square root of salinity
!     zt_degc      : temperature in degrees Celsius
!     zdvi         : volume change for ionization
!     zdki         : compressibility change for ionization
!     zln_kp1_p0   : ln(K_p1) at p_bar = 0
!     zln_kp1_pp   : pressure correction for p_bar /= 0

REAL(KIND=rk) :: zrt, zsqrts, zt_degc, zdvi, zdki
REAL(KIND=rk) :: zlog10_ksparag_p0, zln_ksparag_pp


zsqrts    = SQRT(s)

! log10(Ksp_Arag) for p_bar = 0
zlog10_ksparag_p0 = &
             -171.945_rk - 0.077993_rk*t_k                                     &
            + 2903.293_rk/t_k + 71.595_rk*LOG10(t_k)                            &
            + ( -0.068393_rk + 0.0017276_rk*t_k + 88.135_rk/t_k)*zsqrts        &
            - 0.10018_rk*s + 0.0059415_rk*s*zsqrts


! Pressure correction
zt_degc   = t_k - t_k_zerodegc
zrt       = gasconst_bar_cm3_o_mol_k * t_k

zdvi      =  -48.76_rk + 0.5304_rk*zt_degc  + 2.8_rk
zdki      = (-11.76_rk + 0.3692_rk*zt_degc)*1.0E-03_rk

zln_ksparag_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt


! Final value of Ksp_Arag

ASP_ARAG_MUCC83 = 10._rk**(zlog10_ksparag_p0) * EXP(zln_ksparag_pp)

RETURN

!=======================================================================
 END FUNCTION ASP_ARAG_MUCC83
!=======================================================================







!=======================================================================
FUNCTION A_BTOT_SALIN(s)
!=======================================================================

! Function returns total borate concentration in mol/kg-SW
! given the salinity of a sample

! References: Uppström (1974), cited by  Dickson et al. (2007, chapter 5, p 10)
!             Millero (1982) cited in Millero (1995)
! pH scale  : N/A


IMPLICIT NONE

REAL(KIND=rk) :: A_BTOT_SALIN


! ------------------
! Argument variables
! ------------------

REAL(KIND=rk), INTENT(IN) :: s


A_BTOT_SALIN = 0.000416_rk*(s/35._rk)

RETURN

!=======================================================================
END FUNCTION A_BTOT_SALIN
!=======================================================================



!=======================================================================
FUNCTION A_CATOT_SALIN(s)
!=======================================================================

! Function returns total calcium concentration in mol/kg-SW
! given the salinity of a sample

! References: Culkin (1965)
! pH scale  : N/A


IMPLICIT NONE

REAL(KIND=rk) :: A_CATOT_SALIN


! ------------------
! Argument variables
! ------------------

REAL(KIND=rk), INTENT(IN) :: s


A_CATOT_SALIN = 0.010282_rk*(s/35._rk)


RETURN

!=======================================================================
END FUNCTION A_CATOT_SALIN
!=======================================================================



!=======================================================================
FUNCTION A_FTOT_SALIN(s)
!=======================================================================

! Function returns total calcium concentration in mol/kg-SW
! given the salinity of a sample

! References: Culkin (1965) (???)
! pH scale  : N/A


IMPLICIT NONE

REAL(KIND=rk) :: A_FTOT_SALIN


! ------------------
! Argument variables
! ------------------

REAL(KIND=rk), INTENT(IN) :: s


A_FTOT_SALIN = 0.000068_rk*(s/35._rk)


RETURN

!=======================================================================
END FUNCTION A_FTOT_SALIN
!=======================================================================



!=======================================================================
FUNCTION A_SO4TOT_SALIN(s)
!=======================================================================

! Function returns total sulfate concentration in mol/kg-SW
! given the salinity of a sample

! References: Morris, A.W. and Riley, J.P. (1966) quoted in Handbook (2007)
! pH scale  : N/A


IMPLICIT NONE

REAL(KIND=rk) :: A_SO4TOT_SALIN


! ------------------
! Argument variables
! ------------------

REAL(KIND=rk), INTENT(IN) :: s


!A_SO4TOT_SALIN = 0.028234_rk*(s/35._rk) ! in libthdyct and Thesis
!A_SO4TOT_SALIN = 0.02824_rk*(s/35._rk)                ! Handbook (2007, chap 6, p 10, tab 2, col 3)
A_SO4TOT_SALIN = (0.1400_rk/96.062_rk)*(s/1.80655_rk)  ! Handbook (2007, chap 6, p 10)

RETURN

!=======================================================================
END FUNCTION A_SO4TOT_SALIN
!=======================================================================




!=======================================================================
 FUNCTION ACVT_HSWS_O_HTOT(t_k, s, p_bar)
!=======================================================================

! Function returns the ratio H_SWS/H_Tot as a function of salinity s

! Reference:  Munhoven
! pH scale:   all


IMPLICIT NONE

REAL(KIND=rk) :: ACVT_HSWS_O_HTOT


! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k
REAL(KIND=rk), INTENT(IN) ::  s
REAL(KIND=rk), INTENT(IN) ::  p_bar


! ---------------
! Local variables
! ---------------

!     zso4_tot: total sulfate concentration in mol/kg-SW
!     zf_tot  : total fluoride concentration in mol/kg-SW

REAL(KIND=rk) :: zso4_tot, zf_tot

!-----------------------------------------------------------------------


zso4_tot = A_SO4TOT_SALIN(s)
zf_tot   = A_FTOT_SALIN(s)


ACVT_HSWS_O_HTOT = 1._rk +  (zf_tot*ABETA_HF_DIRI79(t_k, s, p_bar)) &
                           /(1._rk + zso4_tot/AK_HSO4_DICK90(t_k,s, p_bar))

RETURN


!=======================================================================
 END FUNCTION ACVT_HSWS_O_HTOT
!=======================================================================




!=======================================================================
FUNCTION ACVT_HTOT_O_HFREE(t_k, s, p_bar)
!=======================================================================

! Function returns the ratio H_Tot/H_free as a function of salinity s

! Reference:  Munhoven
! pH scale:   N/A


IMPLICIT NONE

REAL(KIND=rk) :: ACVT_HTOT_O_HFREE


! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k
REAL(KIND=rk), INTENT(IN) ::  s
REAL(KIND=rk), INTENT(IN) ::  p_bar


! ---------------
! Local variables
! ---------------

!     zso4_tot: total sulfate concentration in mol/kg-SW

REAL(KIND=rk) :: zso4_tot

!-----------------------------------------------------------------------


zso4_tot = A_SO4TOT_SALIN(s)


ACVT_HTOT_O_HFREE = 1._rk + zso4_tot/AK_HSO4_DICK90(t_k,s, p_bar)

RETURN


!=======================================================================
 END FUNCTION ACVT_HTOT_O_HFREE
!=======================================================================




!=======================================================================
FUNCTION ACVT_HSWS_O_HFREE(t_k, s, p_bar)
!=======================================================================

! Function returns the ratio H_SWS/H_free as a function
! of salinity s

! Reference:  Munhoven
! pH scale:   N/A


IMPLICIT NONE

REAL(KIND=rk) :: ACVT_HSWS_O_HFREE


! ------------------
! Argument variables
! ------------------

!     t_k    : temperature in K
!     s      : salinity
!     p_bar  : applied pressure in bar

REAL(KIND=rk), INTENT(IN) ::  t_k
REAL(KIND=rk), INTENT(IN) ::  s
REAL(KIND=rk), INTENT(IN) ::  p_bar


! ---------------
! Local variables
! ---------------

!     zso4_tot: total sulfate concentration in mol/kg-SW
!     zf_tot  : total fluoride concentration in mol/kg-SW

REAL(KIND=rk) :: zso4_tot, zf_tot

!-----------------------------------------------------------------------


zso4_tot = A_SO4TOT_SALIN(s)
zf_tot   = A_FTOT_SALIN(s)


ACVT_HSWS_O_HFREE = 1._rk + zf_tot*ABETA_HF_DIRI79(t_k, s, p_bar)              &
                          + zso4_tot/AK_HSO4_DICK90(t_k,s, p_bar)

RETURN


!=======================================================================
 END FUNCTION ACVT_HSWS_O_HFREE
!=======================================================================



!=======================================================================
 FUNCTION A_RHOSW1_MUNH97(t_k, s, p_bar)
!=======================================================================

! Function returns first order approximation of \rho in (kg-SW)/(m^3-SW)

! References: Munhoven (1997)
!             after EOS80 (UNESCO, 1981, 1983)

IMPLICIT NONE

REAL(KIND=rk) :: A_RHOSW1_MUNH97

! ------------------
! Argument variables
! ------------------

!     s      : salinity
!     tk     : temperature in K
!     p_bar  : depth in m

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     s0     : 35.5
!     t_k0   : 285.16 K
!     p_bar0 : 300 bar

REAL(KIND=rk), PARAMETER ::     s0 =  35.5_rk
REAL(KIND=rk), PARAMETER ::   t_k0 = 285.16_rk
REAL(KIND=rk), PARAMETER :: p_bar0 = 300.0_rk


A_RHOSW1_MUNH97 = 1039.9044_rk + 0.77629393_rk*(s-s0)                  &
                               - 0.19692738_rk*(t_k-t_k0)                &
                               + 0.044038615_rk*(p_bar-p_bar0)

RETURN

!=======================================================================
 END FUNCTION A_RHOSW1_MUNH97
!=======================================================================



!=======================================================================
 FUNCTION A_RHOSW2_MUNH97(t_k, s, p_bar)
!=======================================================================

! Function returns first order approximation of \rho in (kg-SW)/(m^3-SW)

! References: Munhoven (1997)
!             after EOS80 (UNESCO, 1981, 1983)

IMPLICIT NONE

REAL(KIND=rk) :: A_RHOSW2_MUNH97

! ------------------
! Argument variables
! ------------------

!     s      : salinity
!     tk     : temperature in K
!     p_bar  : depth in m

REAL(KIND=rk), INTENT(IN) ::  t_k, s, p_bar


! ---------------
! Local variables
! ---------------

!     s0     : 35.5
!     t_k0   : 285.16 K
!     p_bar0 : 300 bar

REAL(KIND=rk), PARAMETER ::     s0 =  35.5_rk
REAL(KIND=rk), PARAMETER ::   t_k0 = 285.16_rk
REAL(KIND=rk), PARAMETER :: p_bar0 = 300.0_rk


A_RHOSW2_MUNH97 = 1040.0145_rk                                         &
                   + 0.77629393_rk*(s-s0)                              &
                   - 0.25013591_rk*(t_k-t_k0)                          &
                   + 4.2026266E-02_rk*(p_bar-p_bar0)                   &
                   - 4.7473116E-03_rk*(t_k-t_k0)*(t_k-t_k0)            &
                   - 4.7974224E-06_rk*(p_bar-p_bar0)*(p_bar-p_bar0)    &
                   - 2.1404592E-04_rk*(t_k-t_k0)*(p_bar-p_bar0)

RETURN

!=======================================================================
 END FUNCTION A_RHOSW2_MUNH97
!=======================================================================




!=======================================================================
 SUBROUTINE CHECKCONSTANTS
!=======================================================================

IMPLICIT NONE

! ------------------
! Argument variables
! ------------------

! N/A


! ---------------
! Local variables
! ---------------

!     s      : salinity
!     tk     : temperature in K
!     p_bar  : applied pressure in bar

REAL(KIND=rk) :: t_k, s, p_bar

REAL(KIND=rk) :: zkc0, zkc1, zkc2
REAL(KIND=rk) :: zkb
REAL(KIND=rk) :: zkhf
REAL(KIND=rk) :: zkhso4
REAL(KIND=rk) :: zkp1, zkp2, zkp3
REAL(KIND=rk) :: zksi1
REAL(KIND=rk) :: zkw
REAL(KIND=rk) :: zknh4
REAL(KIND=rk) :: zkh2s

INTEGER, PARAMETER :: logunit = 1


OPEN(logunit,FILE='checkconst.log')

WRITE(logunit,*) 'Checking constant values generated from MOD_CHEMCONST'
WRITE(logunit,*)
WRITE(logunit,*) ' % indicates checking against the Handbook (1994);'
WRITE(logunit,*) ' $ indicates checking against the Lewis and Wallace (1998);'
WRITE(logunit,*) ' * indicates checking against the Handbook (2007);'
WRITE(logunit,*) '   target values are quoted in brackets'
WRITE(logunit,*)
WRITE(logunit,*)
WRITE(logunit,*) ' For S = 35, P = 0 and T/K = 298.15:'
WRITE(logunit,*)


s     = 35._rk
p_bar = 0._rk
t_k   = 298.15_rk


zkc0 = AK_CARB_0_WEIS74(t_k, s)
WRITE(logunit,*) 
WRITE(logunit,*) 'K_0 -- Weiss (1974)'
WRITE(logunit,*) '==================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_0            :', zkc0
WRITE(logunit,*) '   ln(K_0)        :', LOG(zkc0)
WRITE(logunit,*) '   pK_0           :', -LOG10(zkc0)
WRITE(logunit,'("  * ln(K_0)        :", F8.4, " (-3.5617)")') LOG(zkc0)
WRITE(logunit,*) 


zkhso4 = AK_HSO4_DICK90(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_HSO4 -- Dickson (1990) -- pH_free'
WRITE(logunit,*) '==================================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_HSO4         :', zkhso4
WRITE(logunit,*) '   ln(K_HSO4)     :', LOG(zkhso4)
WRITE(logunit,*) '   pK_HSO4        :', -LOG10(zkhso4)
WRITE(logunit,'("  * ln(K_HSO4)     :", F6.2, " (-2.30)")') LOG(zkhso4)
WRITE(logunit,*) 


zkb = AK_BORA_DICK90(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_b -- Dickson (1990) -- pH_tot'
WRITE(logunit,*) '==============================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_b            :', zkb
WRITE(logunit,*) '   ln(K_b)        :', LOG(zkb)
WRITE(logunit,*) '   pK_b           :', -LOG10(zkb)
WRITE(logunit,'("  * ln(K_b)        :", F9.4, " (-19.7964)")') LOG(zkb)
WRITE(logunit,*) 




zkc1 = AK_CARB_1_LUEK00(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_1 -- Luecker et al (2000) -- pH_tot'
WRITE(logunit,*) '====================================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_1            :', zkc1
WRITE(logunit,*) '   ln(K_1)        :', LOG(zkc1)
WRITE(logunit,*) '   pK_1           :', -LOG10(zkc1)
WRITE(logunit,'("  * log10(K_1)     :", F8.4, " (-5.8472)")') LOG10(zkc1)
WRITE(logunit,*) 


zkc2 = AK_CARB_2_LUEK00(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_2 -- Luecker et al (2000) -- pH_tot'
WRITE(logunit,*) '====================================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_2            :', zkc2
WRITE(logunit,*) '   ln(K_2)        :', LOG(zkc2)
WRITE(logunit,*) '   pK_2           :', -LOG10(zkc2)
WRITE(logunit,'("  * log10(K_2)     :", F8.4, " (-8.9660)")') LOG10(zkc2)
WRITE(logunit,*) 




zkc1 = AK_CARB_1_ROYE93(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_1 -- Roy et al (1993) -- pH_tot'
WRITE(logunit,*) '================================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_1            :', zkc1
WRITE(logunit,*) '   ln(K_1)        :', LOG(zkc1)
WRITE(logunit,*) '   pK_1           :', -LOG10(zkc1)
WRITE(logunit,'("  % ln(K_1)        :", F9.4, " (-13.4847)")') LOG(zkc1)
WRITE(logunit,*) 


zkc2 = AK_CARB_2_ROYE93(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_2 -- Roy et al (1993) -- pH_tot'
WRITE(logunit,*) '================================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_2            :', zkc2
WRITE(logunit,*) '   ln(K_2)        :', LOG(zkc2)
WRITE(logunit,*) '   pK_2           :', -LOG10(zkc2)
WRITE(logunit,'("  % ln(K_2)        :", F9.4, " (-20.5504)")') LOG(zkc2)
WRITE(logunit,*) 




zkhf = AK_HF_PEFR87(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_HF -- Perez and Fraga (1987) -- pH_tot'
WRITE(logunit,*) '========================================'
WRITE(logunit,*)
WRITE(logunit,*) '   K_HF           :', zkhf
WRITE(logunit,*) '   ln(K_HF)       :', LOG(zkhf)
WRITE(logunit,*) '   pK_HF          :', -LOG10(zkhf)
WRITE(logunit,'("  * ln(K_HF)       :", F6.2, " (-6.09)")') LOG(zkhf)
WRITE(logunit,*) 


zkp1 = AK_PHOS_1_MILL95(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_P1 -- Millero (1995) -- pH_SWS'
WRITE(logunit,*) '================================'
WRITE(logunit,*)
WRITE(logunit,*) '   K_P1           :', zkp1
WRITE(logunit,*) '   ln(K_P1)       :', LOG(zkp1)
WRITE(logunit,*) '   pK_1           :', -LOG10(zkp1)
WRITE(logunit,'("  * ln(K_P1)-0.015 :", F6.2, " (-3.71)")') LOG(zkp1)-0.015_rk
WRITE(logunit,*) 


zkp2 = AK_PHOS_2_MILL95(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_P2 -- Millero (1995) -- pH_SWS'
WRITE(logunit,*) '================================'
WRITE(logunit,*)
WRITE(logunit,*) '   K_2            :', zkp2
WRITE(logunit,*) '   ln(K_P2)       :', LOG(zkp2)
WRITE(logunit,*) '   pK_2           :', -LOG10(zkp2)
WRITE(logunit,'("  * ln(K_P2)-0.015 :", F8.3, " (-13.727)")') LOG(zkp2)-0.015_rk
WRITE(logunit,*) 


zkp3 = AK_PHOS_3_MILL95(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_P3 -- Millero (1995) -- pH_SWS'
WRITE(logunit,*) '================================'
WRITE(logunit,*)
WRITE(logunit,*) '   K_P3           :', zkp3
WRITE(logunit,*) '   ln(K_P3)       :', LOG(zkp3)
WRITE(logunit,*) '   pK_P3          :', -LOG10(zkp3)
WRITE(logunit,'("  * ln(K_P3)-0.015 :", F7.2, " (-20.24)")') LOG(zkp3)-0.015_rk
WRITE(logunit,*) 


zksi1 = AK_SILI_1_MILL95(t_k, s)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_Si1 -- Millero (1995) -- pH_SWS'
WRITE(logunit,*) '================================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_Si1          :', zksi1
WRITE(logunit,*) '   ln(K_Si1)      :', LOG(zksi1)
WRITE(logunit,*) '   pK_Si1         :', -LOG10(zksi1)
WRITE(logunit,'("  * ln(K_Si1)-0.015:", F7.2, " (-21.61)")') LOG(zksi1)-0.015_rk
WRITE(logunit,*) 


zkw = AK_W_MILL95(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_w -- Millero (1995) -- pH_SWS'
WRITE(logunit,*) '==============================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_w            :', zkw
WRITE(logunit,*) '   ln(K_w)        :', LOG(zkw)
WRITE(logunit,*) '   pK_w           :', -LOG10(zkw)
WRITE(logunit,'("  * ln(K_w)-0.015  :", F8.3, " (-30.434)")') LOG(zkw)-0.015_rk
WRITE(logunit,*) 


zkh2s = AK_H2S_1_MILL95(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_H2S -- Millero (1995) -- pH_SWS'
WRITE(logunit,*) '================================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_H2S          :', zkh2s
WRITE(logunit,*) '   ln(K_H2S)      :', LOG(zkh2s)
WRITE(logunit,*) '   pK_H2S         :', -LOG10(zkh2s)
WRITE(logunit,'("  $ pK_H2S         :", F5.2, " (6.51)")') -LOG10(zkh2s)
WRITE(logunit,*)


zknh4 = AK_AMMO_1_YAMI95(t_k, s, p_bar)

WRITE(logunit,*) 
WRITE(logunit,*) 'K_NH4 -- Yao and Millero (1995) -- pH_SWS'
WRITE(logunit,*) '========================================='
WRITE(logunit,*)
WRITE(logunit,*) '   K_NH4          :', zknh4
WRITE(logunit,*) '   ln(K_NH4)      :', LOG(zknh4)
WRITE(logunit,*) '   pK_NH4         :', -LOG10(zknh4)
WRITE(logunit,'("  $ pK_NH4         :", F5.2, " (9.26)")') -LOG10(zknh4)
WRITE(logunit,*)



CLOSE(logunit)

RETURN

!=======================================================================
 END SUBROUTINE CHECKCONSTANTS
!=======================================================================

END MODULE MOD_CHEMCONST

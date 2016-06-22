!-----------------------------------------------------------------------
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Shamil Yakubov,
!                     Jorn Bruggeman
!-----------------------------------------------------------------------

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:
!
! !INTERFACE:
   module fabm_niva_brom_eqconst
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
    !Gas constant
    real(rk), parameter:: gasconst_bar_cm3_o_mol_k = 83.14472_rk ! Dickson et al. , (2007)
    real(rk), parameter:: t_k_zerodegc = 273.15_rk ! Dickson et al. , (2007)

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_eqconst
!     Variable identifiers
      type (type_diagnostic_variable_id)   :: id_Kc1,id_Kc2,id_Kw,id_Kb,id_Kp1, &
          id_Kp2,id_Kp3,id_Kc0,id_KSi,id_Knh4,id_Kh2s, &
          id_kso4, id_kflu, id_tot_free
      type (type_dependency_id)            :: id_temp,id_salt,id_pres
   contains
      procedure :: initialize
      procedure :: do
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the BROM equilibrium constant model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_eqconst), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s):
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_pres,standard_variables%pressure)

   ! CO2 solubility
   call self%register_diagnostic_variable(self%id_Kc0,'Kc0','-','Henry''s constant')

   ! Carbonic acid
   call self%register_diagnostic_variable(self%id_Kc1,'Kc1','-','[H+][HCO3-]/[H2CO3]')
   call self%register_diagnostic_variable(self%id_Kc2,'Kc2','-','[H+][CO3--]/[HCO3-]')

   ! Water
   call self%register_diagnostic_variable(self%id_Kw, 'Kw','-','[H+][OH-]/H2O')

   ! Boric acid
   call self%register_diagnostic_variable(self%id_Kb, 'Kb','-','[H+][B(OH)4-]/[B(OH)3]')

   ! Phosphoric acid
   call self%register_diagnostic_variable(self%id_Kp1,'Kp1','-','[H+][H2PO4-]/[H3PO4]')
   call self%register_diagnostic_variable(self%id_Kp2,'Kp2','-','[H+][HPO4--]/[H2PO4-]')
   call self%register_diagnostic_variable(self%id_Kp3,'Kp3','-','[H+][PO4---]/[HPO4--]')

   ! Silicic acid
   call self%register_diagnostic_variable(self%id_KSi,'KSi','-','[H+][H3SiO4-]/[Si(OH)4]')

   ! Ammonia
   call self%register_diagnostic_variable(self%id_Knh4,'Knh4','-','[H+][NH3]/[NH4+]')

   ! Hydrogen sulfide
   call self%register_diagnostic_variable(self%id_Kh2s,'Kh2s','-','[H+][HS-]/[H2S]')
  ! call self%register_diagnostic_variable(self%id_Kh2s1,'Kh2s1','-','[H+][HS-]/[H2S]')
  ! call self%register_diagnostic_variable(self%id_Kh2s2,'Kh2s2','-','[H+][S2-]/[HS-]')

   ! Hydrogen sulfate
   call self%register_diagnostic_variable(self%id_kso4,'kso4','-','[H+][HSO4-]/[H2SO4]')

   ! Hydrogen fluoride
   call self%register_diagnostic_variable(self%id_kflu,'kflu','-','[H+][F-]/[HF]')

   ! Ratio H_Tot/H_free
   call self%register_diagnostic_variable(self%id_tot_free,'tot_free','-','ratio H_Tot/H_free')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)

! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_eqconst),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk) :: temp,salt,pres
   real(rk) :: p_bar !applied pressure in bar
   real(rk) :: t_k   !temperature in K
   real(rk) :: Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4, &
       Kh2s, kso4, kflu, tot_free
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_salt,salt)              ! salinity
   _GET_(self%id_pres,pres)              ! pressure in dbar

   p_bar = (pres - 10._rk)/10._rk
   t_k = temp + 273.15_rk

   !calculates equilibrium constants
   !moles per kilogram of solution
   !are expressed in terms of 'total' hydrogen ion concentration
   !tot_free - ratio H_Tot/H_free
   call eq_ph_tot(t_k, salt, p_bar, &
            kc1, kc2, kb, kp1, kp2, kp3, ksi, &
            knh4, kh2s, kso4, kflu, kw, kc0, tot_free)

   ! Transfer all computed equilibrium constants to FABM.
   _SET_DIAGNOSTIC_(self%id_Kc1,  Kc1)
   _SET_DIAGNOSTIC_(self%id_Kc2,  Kc2)
   _SET_DIAGNOSTIC_(self%id_Kw,   Kw)
   _SET_DIAGNOSTIC_(self%id_Kb,   Kb)
   _SET_DIAGNOSTIC_(self%id_Kp1,  Kp1)
   _SET_DIAGNOSTIC_(self%id_Kp2,  Kp2)
   _SET_DIAGNOSTIC_(self%id_Kp3,  Kp3)
   _SET_DIAGNOSTIC_(self%id_Kc0,  Kc0)
   _SET_DIAGNOSTIC_(self%id_KSi,  KSi)
   _SET_DIAGNOSTIC_(self%id_Knh4, Knh4)
   _SET_DIAGNOSTIC_(self%id_Kh2s, Kh2s)
   _SET_DIAGNOSTIC_(self%id_kso4, kso4)
   _SET_DIAGNOSTIC_(self%id_kflu, kflu)
   _SET_DIAGNOSTIC_(self%id_tot_free, tot_free)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!===========================================================================
!
!              eq_ph_tot       !  for pH-TOTAL ! 20.05.2016  Shamil Yakubov
!
!---------------------------------------------------------------------------
    subroutine eq_ph_tot(t_k, s, p_bar, &
            kc1, kc2, kb, kp1, kp2, kp3, ksi, &
            knh4, kh2s, kso4, kflu, kw, kc0, tot_free)

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar
        !dissociation constants mol/kg-SW in total pH scale
        real(rk), intent(out):: kc1 !1st of carbonic acid
        real(rk), intent(out):: kc2 !2nd of carbonic acid
        real(rk), intent(out):: kb  !boric acid
        real(rk), intent(out):: kp1, kp2, kp3 !phosphoric acid
        real(rk), intent(out):: ksi !silicic acid
        real(rk), intent(out):: knh4  !Ammonia
        real(rk), intent(out):: kh2s !Hydrogen sulfide
        real(rk), intent(out):: kso4
        real(rk), intent(out):: kflu
        real(rk), intent(out):: kw  !water
        real(rk), intent(out):: kc0 !Henry's constant
        real(rk), intent(out):: tot_free !ratio H_Tot/H_free

        real(rk):: tot_sws !ratio H_Tot/H_sws

        tot_free = tot_free_ratio(t_k, s, p_bar)
        tot_sws = 1._rk/sws_tot_ratio(t_k, s, p_bar)

        kc1 = carb1_doe94(t_k, s, p_bar)
        kc2 = carb2_doe94(t_k, s, p_bar)
        kb = boric_doe94(t_k, s, p_bar)
        kp1 = phosphoric1_doe94(t_k, s, p_bar)
        kp2 = phosphoric2_doe94(t_k, s, p_bar)
        kp3 = phosphoric3_doe94(t_k, s, p_bar)
        ksi = silicic_doe94(t_k, s)
        knh4 = ammonium_yami95(t_k, s, p_bar)*tot_sws
        kh2s = h2s_mi95(t_k, s, p_bar)*tot_sws
        kso4 = AK_HSO4_DICK90(t_k, s, p_bar)*tot_free
        kflu = flu_pefr87(t_k, s, p_bar)
        kw = water_doe94(t_k, s, p_bar)
        kc0 = carb0_weiss74(t_k, s)

   end subroutine eq_ph_tot

!-----------------------------------------------------------------------

    real(rk) function sws_tot_ratio(t_k, s, p_bar) result(ans)
    !Function returns the ratio H_SWS/H_Tot as a function of salinity s
    !Reference:  Munhoven
    !pH scale:   all
    !t_k    : temperature in K
    !s      : salinity
    !p_bar  : applied pressure in bar

        real(rk), intent(in):: t_k
        real(rk), intent(in):: s
        real(rk), intent(in):: p_bar

        !zso4_tot: total sulfate concentration in mol/kg-SW
        !zf_tot  : total fluoride concentration in mol/kg-SW
        real(rk):: zso4_tot, zf_tot

        !Function returns total sulfate concentration in mol/kg-SW
        zso4_tot = sal_so4(s)
        !Function returns total fluoride concentration in mol/kg-SW
        zf_tot = sal_ftot(s)

        ans = 1._rk +  (zf_tot*ABETA_HF_DIRI79(t_k, s, p_bar)) &
                           /(1._rk + zso4_tot/AK_HSO4_DICK90(t_k,s, p_bar))

    end function sws_tot_ratio

!-----------------------------------------------------------------------

    real(rk) function tot_free_ratio(t_k, s, p_bar) result(ans)
    !Function returns the ratio H_Tot/H_free as a function of salinity s
    !Reference:  Munhoven
    !pH scale:   N/A
    !t_k    : temperature in K
    !s      : salinity
    !p_bar  : applied pressure in bar

        real(rk), intent(in):: t_k
        real(rk), intent(in):: s
        real(rk), intent(in):: p_bar

        !zso4_tot: total sulfate concentration in mol/kg-SW
        real(rk):: zso4_tot

        !Function returns total sulfate concentration in mol/kg-SW
        zso4_tot = sal_so4(s)

        ans = 1._rk + zso4_tot/AK_HSO4_DICK90(t_k,s, p_bar)

    end function tot_free_ratio

!-----------------------------------------------------------------------

    real(rk) function sal_so4(s) result(ans)
    !adopted from Munhoven(2013)
    !Function returns total sulfate concentration in mol/kg-SW
    !given the salinity of a sample
    !References: Morris, A.W. and Riley, J.P. (1966) quoted in Dickson et al. (2007)
    !pH scale  : N/A

        real(rk), intent(in):: s

        ans = (0.1400_rk/96.062_rk)*(s/1.80655_rk)

    end function sal_so4

    real(rk) function sal_ftot(s) result(ans)
    !adopted from Munhoven(2013)
    !Function returns total fluoride concentration in mol/kg-SW
    !given the salinity of a sample
    !References: Culkin (1965) (???)
    !pH scale  : N/A

        real(rk), intent(in):: s

        ans = 0.000068_rk*(s/35._rk)

    end function sal_ftot

!-----------------------------------------------------------------------

    real(rk) function carb1_doe94(t_k, s, p_bar) result(ans)
    !Function calculates first dissociation constant of carbonic acid
    !in mol/kg-sol on the Total pH-scale
    !References: Roy et al. (1993), DOE (1994), Zeebe (2001)
    !            Millero (1995) pressure correction
    !pH scale:   Total

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_kc1_p0      : ln_kc1 at p_bar = 0
        !ln_kc1_pp      : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_kc1_p0, ln_kc1_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        real(rk):: TKR, TLOG, SQ, SLOG!, pk

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk
        !p_bar = 300._rk

        TKR  = 1./t_k
        TLOG = LOG(t_k)
        SQ   = SQRT(s)
        SLOG = LOG(1.0e0_rk-0.001005e0_rk*s)

        ln_kc1_p0 = -2307.1266_rk*TKR + 2.83655_rk - 1.5529413_rk*TLOG &
              + (-4.0484_rk*TKR - 0.20760841_rk - 0.00654208_rk*s)*SQ  &
              + 0.08468345_rk*s + SLOG

        !Pressure correction
        zt_degc     = t_k - t_k_zerodegc
        zrt         = gasconst_bar_cm3_o_mol_k * t_k
        zdvi        = -25.50_rk + zt_degc*(0.1271_rk + zt_degc*0.0E-03_rk)
        zdki        = (-3.08_rk + zt_degc*0.0877_rk)*1.0E-03_rk
        ln_kc1_pp   = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_kc1_p0 + ln_kc1_pp)

        !for T=25C, S=35psu, P=0bars:    pk = 5.8563
        !for T=25C, S=35psu, P=300bars:  pk = 5.7397
        !CHECKED 13_05_2016 Shamil
        !pk = -LOG10(ans)

    end function carb1_doe94

!-----------------------------------------------------------------------

    real(rk) function carb2_doe94(t_k, s, p_bar) result(ans)
    !Function calculates second dissociation constant of carbonic acid
    !in mol/kg-sol on the Total pH-scale
    !References: Roy et al. (1993), DOE (1994), Zeebe (2001)
    !            Millero (1995) pressure correction
    !pH scale:   Total

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_kc2_p0      : ln_kc2 at p_bar = 0
        !ln_kc2_pp      : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_kc2_p0, ln_kc2_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        real(rk):: TKR, TLOG, SQ, SLOG!, pk

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk
        !p_bar = 300._rk

        TKR  = 1./t_k
        TLOG = LOG(t_k)
        SQ   = SQRT(s)
        SLOG = LOG(1.0e0_rk-0.001005e0_rk*s)

        ln_kc2_p0 = -3351.6106_rk*TKR - 9.226508_rk - 0.2005743_rk*TLOG &
              + (-23.9722_rk*TKR - 0.106901773_rk - 0.00846934_rk*s)*SQ &
              + 0.1130822_rk*s + SLOG

        !Pressure correction
        zt_degc     = t_k - t_k_zerodegc
        zrt         = gasconst_bar_cm3_o_mol_k * t_k
        zdvi        = -15.82_rk + zt_degc*(-0.0219_rk + zt_degc*0.0E-03_rk)
        zdki        = (1.13_rk + zt_degc*(-0.1475_rk))*1.0E-03_rk
        ln_kc2_pp   = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_kc2_p0 + ln_kc2_pp)

        !for T=25C, S=35psu, P=0bars:    pk = 8.9249
        !for T=25C, S=35psu, P=300bars:  pk = 8.8409
        !CHECKED 13_05_2016 Shamil
        !pk = -LOG10(ans)

    end function carb2_doe94

!-----------------------------------------------------------------------

    real(rk) function boric_doe94(t_k, s, p_bar) result(ans)
    !Function calculates dissociation constant of boric acid
    !in mol/kg-sol on the Total pH-scale
    !References: Dickson (1990), DOE (1994), Zeebe (2001)
    !            Millero (1995) pressure correction
    !pH scale:   Total

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_kb_p0       : ln_kb at p_bar = 0
        !ln_kb_pp       : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_kb_p0, ln_kb_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        real(rk):: TKR, TLOG, SQ!, pk

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk
        !p_bar = 300._rk

        TKR  = 1./t_k
        TLOG = LOG(t_k)
        SQ   = SQRT(s)

        ln_kb_p0 = (-8966.90_rk + (1.728_rk*s - 2890.53_rk)*SQ     &
                          - (77.942_rk + 0.0996_rk*s)*s)*TKR       &
              + (148.0248_rk + 137.1942_rk*SQ + 1.62142_rk*s)      &
              - ( 24.4344_rk +  25.0850_rk*SQ + 0.24740_rk*s)*TLOG &
              + 0.053105_rk*SQ*t_k

        !Pressure correction
        zt_degc     = t_k - t_k_zerodegc
        zrt         = gasconst_bar_cm3_o_mol_k * t_k
        zdvi        = -29.48_rk + zt_degc*(0.1622_rk + zt_degc*2.608E-03_rk)
        zdki        = (-2.84_rk + zt_degc*(0._rk))*1.0E-03_rk
        ln_kb_pp   = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_kb_p0 + ln_kb_pp)

        !for T=25C, S=35psu, P=0bars:    pk = 8.5975
        !for T=25C, S=35psu, P=300bars:  pk = 8.4746
        !CHECKED 13_05_2016 Shamil
        !pk = -LOG10(ans)

    end function boric_doe94

!-----------------------------------------------------------------------

    real(rk) function phosphoric1_doe94(t_k, s, p_bar) result(ans)
    !Function calculates first dissociation constant of phosphoric acid
    !in mol/kg-sol on the Total pH-scale
    !References: DOE (1994), Zeebe (2001)
    !            Millero (1995) pressure correction
    !pH scale:   Total

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_kp1_p0      : ln_kp1 at p_bar = 0
        !ln_kp1_pp      : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_kp1_p0, ln_kp1_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        real(rk):: TKR, TLOG, SQ!, pk

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk
        !p_bar = 300._rk

        TKR  = 1./t_k
        TLOG = LOG(t_k)
        SQ   = SQRT(s)

        ln_kp1_p0 = -4576.752_rk*TKR + 115.525_rk - 18.453_rk*TLOG &
               + (0.69171_rk - 106.736_rk*TKR)*SQ &
               - (0.01844_rk + 0.65643_rk*TKR)*s

        !Pressure correction
        zt_degc     = t_k - t_k_zerodegc
        zrt         = gasconst_bar_cm3_o_mol_k * t_k
        zdvi        = -14.51_rk + zt_degc*(0.1211_rk + zt_degc*(-0.321E-03_rk))
        zdki        = (-2.67_rk + zt_degc*0.0427_rk)*1.0E-03_rk
        ln_kp1_pp   = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_kp1_p0 + ln_kp1_pp)

        !for T=25C, S=35psu, P=0bars:    pk = 1.61
        !for T=25C, S=35psu, P=300bars:  pk = unknown (1.55 in here)
        !CHECKED 13_05_2016 Shamil
        !pk = -LOG10(ans)

    end function phosphoric1_doe94

!-----------------------------------------------------------------------

    real(rk) function phosphoric2_doe94(t_k, s, p_bar) result(ans)
    !Function calculates second dissociation constant of phosphoric acid
    !in mol/kg-sol on the Total pH-scale
    !References: DOE (1994), Zeebe (2001)
    !            Millero (1995) pressure correction
    !pH scale:   Total

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_kp2_p0      : ln_kp2 at p_bar = 0
        !ln_kp2_pp      : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_kp2_p0, ln_kp2_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        real(rk):: TKR, TLOG, SQ!, pk

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk
        !p_bar = 300._rk

        TKR  = 1./t_k
        TLOG = LOG(t_k)
        SQ   = SQRT(s)

        ln_kp2_p0 = -8814.715_rk*TKR + 172.0883_rk - 27.927_rk*TLOG &
               + (1.35660_rk - 160.340_rk*TKR)*SQ &
               - (0.05778_rk - 0.37335_rk*TKR)*s

        !Pressure correction
        zt_degc     = t_k - t_k_zerodegc
        zrt         = gasconst_bar_cm3_o_mol_k * t_k
        zdvi        = -23.12_rk + zt_degc*(0.1758_rk + zt_degc*(-2.647E-03_rk))
        zdki        = (-5.15_rk + zt_degc*0.09_rk)*1.0E-03_rk
        ln_kp2_pp   = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_kp2_p0 + ln_kp2_pp)

        !for T=25C, S=35psu, P=0bars:    pk = 5.96
        !for T=25C, S=35psu, P=300bars:  pk = unknown (5.86 in here)
        !CHECKED 13_05_2016 Shamil
        !pk = -LOG10(ans)

    end function phosphoric2_doe94

!-----------------------------------------------------------------------

    real(rk) function phosphoric3_doe94(t_k, s, p_bar) result(ans)
    !Function calculates third dissociation constant of phosphoric acid
    !in mol/kg-sol on the Total pH-scale
    !References: DOE (1994), Zeebe (2001)
    !            Millero (1995) pressure correction
    !pH scale:   Total

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_kp3_p0      : ln_kp3 at p_bar = 0
        !ln_kp3_pp      : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_kp3_p0, ln_kp3_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        real(rk):: TKR, TLOG, SQ!, pk

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk
        !p_bar = 300._rk

        TKR  = 1./t_k
        TLOG = LOG(t_k)
        SQ   = SQRT(s)

        ln_kp3_p0 = -3070.75_rk*TKR - 18.141_rk &
               + (2.81197_rk + 17.27039_rk*TKR)*SQ &
               - (0.09984_rk + 44.99486_rk*TKR)*s

        !Pressure correction
        zt_degc     = t_k - t_k_zerodegc
        zrt         = gasconst_bar_cm3_o_mol_k * t_k
        zdvi        = -26.57_rk + zt_degc*(0.2020_rk + zt_degc*(-3.042E-03_rk))
        zdki        = (-4.08_rk + zt_degc*0.0714_rk)*1.0E-03_rk
        ln_kp3_pp   = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_kp3_p0 + ln_kp3_pp)

        !for T=25C, S=35psu, P=0bars:    pk = 8.79
        !for T=25C, S=35psu, P=300bars:  pk = unknown (8.67 in here)
        !CHECKED 13_05_2016 Shamil
        !pk = -LOG10(ans)

    end function phosphoric3_doe94

!-----------------------------------------------------------------------

    real(rk) function silicic_doe94(t_k, s) result(ans)
    !Function calculates dissociation constant of silicic acid
    !in mol/kg-sol on the Total pH-scale
    !References: DOE (1994), Zeebe (2001)
    !pH scale:   Total

        !t_k    : temperature in K
        !s      : salinity
        real(rk), intent(in):: t_k, s

        !ln_ksi_p0      : ln_ksi at p_bar = 0
        real(rk):: ln_ksi_p0

        real(rk):: TLOG, SLOG, I!, pk

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk

        TLOG = LOG(t_k)
        SLOG = LOG(1.0e0_rk-0.001005e0_rk*s)
        I    = 19.924_rk*s/(1000._rk-1.005_rk*s)

        ln_ksi_p0 = -8904.2_rk/t_k+117.385_rk-19.334_rk*TLOG &
           + (-458.79_rk/t_k+3.5913_rk)*SQRT(I) &
           + (188.74_rk/t_k-1.5998_rk)*I+(-12.1652_rk/t_k+0.07871_rk)*I*I + SLOG

        ans = EXP(ln_ksi_p0)

        !for T=25C, S=35psu, P=0bars:    pk = 9.38
        !CHECKED 13_05_2016 Shamil
        !pk = -LOG10(ans)

    end function silicic_doe94

!-----------------------------------------------------------------------

    real(rk) function ammonium_yami95(t_k, s, p_bar) result(ans)
    !adopted from Munhoven(2013)
    !Function calculates dissociation constant of ammonium
    !in mol/kg-sol on the Total pH-scale
    !References: Yao and Millero (1995)
    !            Millero (1995) pressure correction
    !pH scale:   SWS

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_knh4_p0     : ln_knh4 at p_bar = 0
        !ln_knh4_pp     : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_knh4_p0, ln_knh4_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        ln_knh4_p0 = -0.25444_rk -  6285.33_rk/t_k + 0.0001635_rk*t_k &
               + ( 0.46532_rk - 123.7184_rk/t_k) * SQRT(s)            &
               + (-0.01992_rk +  3.17556_rk/t_k) * s

        !Pressure correction
        zt_degc     = t_k - t_k_zerodegc
        zrt         = gasconst_bar_cm3_o_mol_k * t_k
        zdvi        =  -26.43_rk + zt_degc*(0.0889_rk - zt_degc*0.905E-03_rk)
        zdki        = ( -5.03_rk + zt_degc*0.0814_rk)*1.0E-03_rk
        ln_knh4_pp  = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_knh4_p0 + ln_knh4_pp)

    end function ammonium_yami95

!-----------------------------------------------------------------------

    real(rk) function h2s_mi95(t_k, s, p_bar) result(ans)
    !adopted from Munhoven(2013)
    !Function calculates dissociation constant of hydrogen sulfide
    !in mol/kg-sol
    !References: Millero et al. (1988) (cited by Millero (1995)
    !            Millero (1995) for pressure correction
    !pH scale  : - SWS (according to Yao and Millero, 1995,
    !              p. 82: "refitted if necessary")
    !            - Total (according to Lewis and Wallace, 1998)
    !Note      : we stick to SWS here for the time being
    !Note      : the fits from Millero (1995) and Yao and Millero (1995)
    !            derive from Millero et al. (1998), with all the coefficients
    !            multiplied by -ln(10

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_kh2s_p0     : ln_kh2s at p_bar = 0
        !ln_kh2s_pp     : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_kh2s_p0, ln_kh2s_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        ln_kh2s_p0 = 225.838_rk-13275.3_rk/t_k &
                   - 34.6435_rk*LOG(t_k)       &
                   + 0.3449_rk*SQRT(s)-0.0274_rk*s

        !Pressure correction
        zt_degc     = t_k - t_k_zerodegc
        zrt         = gasconst_bar_cm3_o_mol_k * t_k
        zdvi        =  -14.80_rk + zt_degc*(0.0020_rk - zt_degc*0.400E-03_rk)
        zdki        = (  2.89_rk + zt_degc*0.054_rk)*1.0E-03_rk
        ln_kh2s_pp  = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_kh2s_p0 + ln_kh2s_pp)

    end function h2s_mi95

!-----------------------------------------------------------------------

    real(rk) FUNCTION ABETA_HF_DIRI79(t_k, s, p_bar) result(ans)
    !adopted from Munhoven(2013)
    !Function calculates association constant \beta_{HF} [(mol/kg-SW)^{-1}]
    !in (mol/kg-SW)^{-1}, where
    !  \beta_{HF} = \frac{ [HF] }{ [H^{+}] [F^{-}] }
    !References: Dickson and Riley (1979)
    !            Millero (1995) for pressure correction
    !pH scale  : free
    !Note      : converted here from mol/kg-H2O to mol/kg-SW

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        !zionst         : ionic strength in mol/-kg-H2O
        !zcvt_to_kgsw   : mass of pure water in 1kg of seawater as a fct. of salinity
        !zln_bhf_p0     : \beta_HF at p_bar = 0
        !zln_khf_pp     : pressure correction for k_HF = 1/\beta_HF at p_bar /= 0

        real(rk):: zrt, zt_degc, zdvi, zdki
        real(rk):: zionst, zcvt_to_kgsw
        real(rk):: zln_bhf_p0, zln_khf_pp
        !ionic strength in mol/kg-SW, for given salinity
        real(rk):: ion_strength

        !\beta_HF at p_bar = 0
        zcvt_to_kgsw = 1._rk - 0.001005_rk*s !Dickson et al.  (2007)
        ion_strength = (0.019924D+00*s) !Dickson et al.  (2007)
        zionst       = ion_strength/zcvt_to_kgsw
        zln_bhf_p0   = -1590.2_rk/t_k + 12.641_rk - 1.525_rk*SQRT(zionst)

        ! Pressure correction
        zt_degc      = t_k - t_k_zerodegc
        zrt          = gasconst_bar_cm3_o_mol_k * t_k
        zdvi         =   -9.78_rk + zt_degc*(-0.0090_rk - zt_degc*0.942E-03_rk)
        zdki         = ( -3.91_rk + zt_degc*0.054_rk)*1.0E-03_rk
        zln_khf_pp   = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ! Final \beta_HF value
        ! notice that  ln(k_HF(P)) = ln(k_HF(0)) + zln_khf_pp
        !         <=>  -ln(\beta_HF(P)) = -ln(\beta_HF(0)) + zln_khf_pp
        !         <=>   ln(\beta_HF(P)) =  ln(\beta_HF(0)) - zln_khf_pp
        ans = EXP(zln_bhf_p0 - zln_khf_pp ) / zcvt_to_kgsw

    END FUNCTION ABETA_HF_DIRI79

!-----------------------------------------------------------------------

    real(rk) FUNCTION AK_HSO4_DICK90(t_k, s, p_bar) result(ans)
    !adopted from Munhoven(2013)
    !Function returns the dissociation constant of hydrogen sulfate (bisulfate)
    !References: Dickson (1990) -- also Dickson et al.  (2007)
    !            Millero (1995) for pressure correction
    !pH scale  : free
    !Note      : converted here from mol/kg-H2O to mol/kg-SW

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        !zionst         : ionic strength in mol/-kg-H2O
        !zsqrti         : square root og ion strength
        !zcvt_to_kgsw   : mass of pure water in 1kg of seawater as a fct. of salinity
        !zln_khso4_p0   : K_HSO4 at p_bar = 0
        !zln_khso4_pp   : pressure correction for p_bar /= 0
        real(rk):: zrt, zt_degc, zdvi, zdki
        real(rk):: zcvt_to_kgsw, zionst, zsqrti
        real(rk):: zln_khso4_p0, zln_khso4_pp
        !ionic strength in mol/kg-SW, for given salinity
        real(rk):: ion_strength!, ln_khso4

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk
        !p_bar = 0._rk

        zcvt_to_kgsw = 1._rk - 0.001005_rk*s !Dickson et al.  (2007)
        ion_strength = (0.019924D+00*s) !Dickson et al.  (2007)
        zionst       = ion_strength/zcvt_to_kgsw
        zsqrti       = SQRT(zionst)
        zln_khso4_p0 =     -4276.1_rk/t_k + 141.328_rk -  23.093_rk*LOG(t_k)&
        + (-13856._rk/t_k  +  324.57_rk -  47.986_rk*LOG(t_k)) * zsqrti &
        + ( 35474._rk/t_k  -  771.54_rk + 114.723_rk*LOG(t_k)) * zionst &
        - (  2698._rk/t_k)*zsqrti * zionst                              &
        + (  1776._rk/t_k)*zionst * zionst

        !Pressure correction
        zt_degc      = t_k - t_k_zerodegc
        zrt          = gasconst_bar_cm3_o_mol_k * t_k
        zdvi         =  -18.03_rk + zt_degc*(0.0466_rk + zt_degc*0.316E-03_rk)
        zdki         = (-4.53_rk + zt_degc*0.0900_rk)*1.0E-03_rk
        zln_khso4_pp = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = zcvt_to_kgsw * EXP(zln_khso4_p0 + zln_khso4_pp)

        !for T=25C, S=35psu, P=0bars: ln_khso4 = -2.3 Dickson et al.  2007
        !CHECKED 16_05_2016 Shamil
        !ln_khso4 = zln_khso4_p0 + log(zcvt_to_kgsw)

    END FUNCTION AK_HSO4_DICK90

    !-----------------------------------------------------------------------

    real(rk) function flu_pefr87(t_k, s, p_bar) result(ans)
    !adopted from Munhoven(2013)
    !Function calculates dissociation constant of hydrogen fluoride
    !in mol/kg-sw
    !References: Perez and Fraga (1987)
    !            Millero (1995) for pressure correction
    !pH scale  : Total (according to Dickson et al. , 2007)

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_khf_p0      : ln_khf at p_bar = 0
        !ln_khf_pp      : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_khf_p0, ln_khf_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk
        !p_bar = 0._rk

        ln_khf_p0 = 874._rk/t_k-9.68_rk+0.111_rk*SQRT(s)

        !Pressure correction
        zt_degc    = t_k - t_k_zerodegc
        zrt        = gasconst_bar_cm3_o_mol_k * t_k
        zdvi       =   -9.78_rk + zt_degc*(-0.0090_rk - zt_degc*0.942E-03_rk)
        zdki       = ( -3.91_rk + zt_degc*0.054_rk)*1.0E-03_rk
        ln_khf_pp  = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_khf_p0 + ln_khf_pp)

        !for T=25C, S=35psu, P=0bars: ln_khf_p0 = -6.09 (Dickson et al., 2007)
        !CHECKED 16_05_2016 Shamil

    end function flu_pefr87

!-----------------------------------------------------------------------

    real(rk) function water_doe94(t_k, s, p_bar) result(ans)
    !Function calculates water dissociation constant Kw
    !in (mol/kg-sol)^2
    !References: DOE (1994), Zeebe (2001)
    !            Millero (1995) for pressure correction
    !pH scale  : Total

        !t_k    : temperature in K
        !s      : salinity
        !p_bar  : applied pressure in bar
        real(rk), intent(in):: t_k, s, p_bar

        !ln_kw_p0       : ln_kw at p_bar = 0
        !ln_kw_pp       : pressure correction for p_bar /= 0
        !zrt            : R*t_k, R in bar*cm3/(mol*K)
        !zt_degc        : temperature in degrees Celsius
        !zdvi           : volume change for ionization
        !zdki           : compressibility change for ionization
        real(rk):: ln_kw_p0, ln_kw_pp
        real(rk):: zt_degc, zrt, zdvi, zdki

        real(rk):: TKR, TLOG, SQ!, pk

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk
        !p_bar = 300._rk

        TKR  = 1./t_k
        TLOG = LOG(t_k)
        SQ   = SQRT(s)

        ln_kw_p0 = -13847.26_rk*TKR + 148.96502_rk - 23.6521_rk*TLOG &
              + (118.67_rk*TKR - 5.977_rk + 1.0495_rk*TLOG)*SQ - 0.01615_rk*s

        !Pressure correction
        zt_degc   = t_k - t_k_zerodegc
        zrt       = gasconst_bar_cm3_o_mol_k * t_k
        zdvi      =   -25.6_rk + zt_degc*(0.2324_rk - zt_degc*3.6246E-03_rk)
        zdki      = ( -5.13_rk + zt_degc*0.0794_rk)*1.0E-03_rk
        ln_kw_pp  = (-zdvi + zdki*p_bar/2._rk)*p_bar/zrt

        ans = EXP(ln_kw_p0 + ln_kw_pp)

        !for T=25C, S=35psu, P=0bars:    pk = 13.2173
        !for T=25C, S=35psu, P=300bars:  pk = 13.1039
        !CHECKED 16_05_2016 Shamil
        !pk = -LOG10(ans)

    end function water_doe94

!-----------------------------------------------------------------------

    real(rk) function carb0_weiss74(t_k, s) result(ans)
    !Function calculates K0 in (mol/kg-SW)/atmosphere
    !References: Weiss (1979), DOE (1994), Zeebe (2001)
    !pH scale  : N/A

        !t_k    : temperature in K
        !s      : salinity
        real(rk), intent(in):: t_k, s

        !ln_kc0_p0       : ln_kc0 at p_bar = 0
        real(rk):: ln_kc0_p0!, pk

        !t_k = 25._rk + 273.15_rk
        !s = 35._rk

        ln_kc0_p0 = -60.2409_rk+9345.17_rk/t_k+23.3585_rk*log((t_k)/100._rk) &
            +s*(0.023517_rk-0.023656_rk*t_k/100._rk &
            +0.0047036_rk*((t_k/100._rk)*(t_k/100._rk)))

        ans = EXP(ln_kc0_p0)

        !for T=25C, S=35psu, P=0bars:    pk = 1.5468
        !CHECKED 16_05_2016 Shamil
        !pk = -LOG10(ans)

    end function carb0_weiss74

end module

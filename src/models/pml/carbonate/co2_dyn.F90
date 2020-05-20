#include "fabm_driver.h"

! This file is based upon MEECE deliverable D2.2, see http://www.meece.eu.
!
! Minor changes were made to make it fit within FABM:
!
! - inserted a reference to the FABM include file (fabm_driver.h) to permit use of FABM preprocessor
!   macros (see below)
! - removed the program - end program section to allow reuse as software library in other programs.
! - relaxed lower bound on temperature, from 0 to -4 degrees Celsius, to permit use in cold waters
!   [after consultation with Jerry Blackford]
! - replaced fixed double precision type (real*8) by preprocessor macro real(rk).
! - replaced stop statements by calls to fatal_error.
!
! Jorn Bruggeman, 4 March 2011

! This file contains a set of FORTRAN subroutines that calculate the carbonate system
! at any given point in marine space time, given values for
! temperature, salinity, DIC, depth (pressure).
! This is essentially an implimentation of the Haltafall speciation code
! (Ingri et al 1967, Talanta 14, 1261 - if it ain't broke don't fix it)
! Another routine calulates the air sea exchange of CO2 given wind speed and atmospheric pCO2.
! Code developed by Jerry blackford and others at PML, based on pre-existing code.
! We accept no liability for errors or inaccuracies.
! See Zeebe & Wolf-Gladrow, 2001. CO2 in seawater: equilibrium, kinetics and isotopes.
! Elsevier Oceanography Series 65, 346. for a reasonable overview.
! Many other packages exist, replicating the same functionality in different languages.
! See http://cdiac.ornl.gov/oceans/co2rprt.html (CO2sys)
! or  http://neon.otago.ac.nz/research/mfc/people/keith_hunter/software/swco2/
! reference for prior usage of this code: Blackford & Gilbert, 2007. J Mar Sys 64, 229-241.

! Modifications
! 17/02/2010. Added conversion factor from per m3 to per kg (line 108-133)
! 17/02/2010. Update calculation of K1, K2, Kb to make consistant with the OCMIP protocols.

   subroutine CO2_dynamics(T,S,Z,DIC,pco2w,TA,ph,carba,bicarb,carb,henry,om_cal,om_arg,TCO2,dcf)

   use fabm_types,only:rk

   IMPLICIT NONE
! INPUT PARAMETERS:
      real(rk)   :: T, S, Z, DIC
! T: temperature (C)
! S: salinity (psu)
! Z: depth (metres)
! DIC: total dissolved inorganic carbon (also commonly called TC total carbon), (mmol.m-3)
! pCO2a: partial pressure of co2 in atmosphere (eg ~390 in 2010, 280=preindustrial)

! LOCAL VARIABLES:
      real(rk)   :: TCO2, TA, ph, pco2w
      real(rk)   :: pco2water,fairco2, henry
      real(rk)   :: ca, bc, cb, carba, bicarb, carb, om_cal, om_arg
      real(rk)   :: a,b,c,dcf

!!!!!!!! provide an expression for total alkalinity.!!!!!!!!!!!!!!
! Note, for oceanic regimes there is generally a well constrained relationship between salinity and
! total alkalinity as TA is conservative. However the slope and intercept vary according to region.
! For coastal regions, such relationships are unsatisfactory.
! If you have the data then one possibility is to introduce riverine alkalinity as a seperate tracer.
! note that uptake/release of nutrients and Calcium carbonate also make small modifications to local TA.
! See Zeebe & Wolf-Gladrow; Millero et al 1998, Mar Chem 60, 111-130; Lee et al, 2006, GRL 33, L19605,
! Wolf-Gladrow et al, Mar Chem 106 (2007).

! deriving alkalinity (umol/kg) from salinity
!          TA = 66.96*S - 36.803 ! Nordic, N Atlantic, Bellerby et al 2005.
!          TA = 520.1 + 51.24*S  ! Atlantic (Millero 1998)
!          TA = 399.0 + 54.629*S ! Pacific (Millero, 1998)
!          TA = -114.0 + 68.80*S ! Indian (Millero, 1998)
          TA = 520.1 + 51.24*S  ! Atlantic (Millero 1998)

! Adjust to correct units.
! Haltafall uses mol/kg rather than umol/kg necessitating a scaling factor of /1.0D6
! DIC (mmol/m3) needs to be converted to umol/kg via the calculation of water density at prevailing T&S
! sea-water density (Millero & Poisson, Deep-Sea Research, 1981, also known as UNESCO, 1981)
! with T: Temperature in degree Celsius; S: Salinity in practical units, density in kg/m3
! valid for 0<T<40 and 0.5<S<43. Density Conversion factor (dcf) = * density * 1.0D3

          a=  8.24493d-1 - 4.0899d-3*T +  7.6438d-5*T**2 - 8.2467d-7*T**3 + 5.3875d-9*T**4
          b= -5.72466d-3 + 1.0227d-4*T - 1.6546d-6*T**2
          c= 4.8314d-4
          dcf= (999.842594 + 6.793952d-2*T- 9.095290d-3*T**2 + 1.001685d-4*T**3 &
                      - 1.120083d-6*T**4 + 6.536332d-9*T**5+a*S+b*S**1.5+c*S**2)/1.0D3

          TA = TA / (1.0D6)
          TCO2  = DIC / (1.0D6*dcf)

!..call the parent routine for the carbonate system
          CALL CO2DYN ( TCO2, TA, T, S, PCO2WATER, pH, HENRY, ca, bc, cb)

! Adjust outputs back to units used in the parent model code (e.g. mmol/m3) if appropriate

          PCO2W = PCO2WATER*1.0D6   ! partial pressure of co2 in water
          TA = TA*(1.0D6)           ! total alkalinity (umol/kg)
          CarbA= ca*(1.0D6*dcf)     ! carbonic acid concentration (mmol/m3)
          Bicarb=bc*(1.0D6*dcf)     ! bicarbonate ion concentration (mmol/m3)
          Carb = cb*(1.0D6*dcf)     ! carbonate ion concentration (mmol/m3)
          TCO2 = TCO2*1.0D6         ! total C or DIC in units of umol/kg

!Call carbonate saturation state subroutine to calculate calcite and aragonite calcification states

          CALL CaCO3_Saturation ( T, S, Z, cb, Om_cal, Om_arg)

! outputs are
!         Om_cal        calcite saturation state
!         Om_arg        aragonite saturation state

   RETURN
   END SUBROUTINE CO2_dynamics

   SUBROUTINE Air_sea_exchange (T, Wnd, pCO2w, pCO2a, Henry, dcf, flux)
!  this routine should be called for the surface box only
!  Uses the Nightingale and Liss parameterisation (GBC 14, 373-388 (2000).
!  SC is the Schmidt number from Wanninkhof (1992, JGR 97,7373-7382), the Nightingale et al (2000)
!  transfer velocity being for Sc=600 and the dependence of tranfer velocity on Sc coming
!  from Jahne et al (1987, J. Geophys Res 92, 1937-1949).

!  Inputs
!  pCO2w    partial pressure of CO2 in the water (from carbonate system subroutine call)
!  pCO2a    partial pressure of CO2 in the atmosphere (usually external forcing).
!  T        temperature (C)
!  Wnd      wind speed, metres
!  Henry    henry's constant
!  density  the density of water for conversion between mmol/m3 and umol/kg

!  Outputs are
!  flux     flux of CO2 in mmol C /m2/d
!           +ve is in-gassing (air to sea), -ve is outgassing (sea to air).
   use fabm_types

   IMPLICIT NONE
      real(rk)   :: T, wnd, pco2w, pco2a, henry, dcf    ! INPUT PARAMETERS:
      real(rk)   :: sc, fwind    ! LOCAL VARIABLES:
      real(rk)   :: flux    ! OUTPUT Variables

! calculate the scmidt number and unit conversions
          sc=2073.1-125.62*T+3.6276*T**2.0-0.0432190*T**3.0
          fwind =  (0.222_rk * wnd**2_rk + 0.333_rk * wnd)*(sc/660._rk)**(-0.5_rk)
          fwind=fwind*24._rk/100._rk   ! convert to m/day

! flux depends on the difference in partial pressures, wind and henry
! here it is rescaled to mmol/m2/d
          flux = fwind * HENRY * ( PCO2A - PCO2W ) * dcf

  RETURN
  END SUBROUTINE Air_sea_exchange

      SUBROUTINE CO2dyn ( TCO2, TA, T, S, pco2, ph, henry, ca, bc, cb)

!.......................................................................
!     This subroutine acts as an interface to the Haltafall iteration, setting options etc.
!.......................................................................
      use fabm_types

      IMPLICIT NONE

      real(rk) PRSS, PH, AKVAL, CONCS,                       &
     &                 TCO2, TA, T, S, PCO2,            &
     &                 SOLBTY, CCO2,                                 &
     &                 A1, A2, A3, B1, B2, B3, TK, TK1, SOL1, SOL2,  &
     &                 HENRY, ca, bc, cb
      INTEGER MCONC, MKVAL, ICONST, ICALC

! INPUT PARAMETERS:
      PARAMETER ( MCONC = 9,MKVAL = 4 )

      DIMENSION AKVAL(MKVAL), CONCS(MCONC)

      ICONST = 6
      PRSS = 1.0_rk
      CONCS(1) = TCO2
      CONCS(2) = TA
      ICALC = 1

      CALL POLYCO(PRSS,T,S,CONCS,MCONC,AKVAL,MKVAL,ICALC,ICONST)

      PCO2   = CONCS(3)
      PH = CONCS(4)
      ca = CONCS(5)
      bc = CONCS(6)
      cb = CONCS(7)
      HENRY = AKVAL(1)

      RETURN
      END SUBROUTINE

      SUBROUTINE POLYCO(PD,TD,SD,CONCS,NCONC,AKVAL,NKVAL,ICALC,ICONST)
!     -----------------------------------------------------------------
! MASTER SUBROUTINE FOR CALCULATION OF THE CO2 SYSTEM THERMODYNAMICS
!

! EXPLANATION OF POLYCO PARAMETERS
!       P - PRESSURE IN ATMOSPHERES (P<>1 NOT YET CODED)
!       T - TEMPERATURE IN DEG.C
!       S - SALINITY IN PPT
!      CONCS(1) - TOTAL C (MOL/KG)
!      CONCS(2) - TOTAL ALKALINITY (MOL/KG)
!      CONCS(3) - PCO2 (ATM)
!      CONCS(4) - PH
!      CONCS(5) - {H2CO3} (MOL/KG)
!      CONCS(6) - {HCO3} (MOL/KG)
!      CONCS(7) - {CO3} (MOL/KG)
!        CONCS(8) - CARBONATE ALKALINITY  ) FOR ICONST = 4,5,6
!        CONCS(9) - BORATE ALKALINITY     )       ONLY
!         NCONC - SIZE OF CONCS ARRAY (7 FOR ICONST=1,2,3; 9 FOR ICONST
!     AKVAL(1) - KP (HENRY'S LAW CONSTANT) (MOL/KG/ATM)
!     AKVAL(2) - K1C (H2CO3 DISSOCIATION) (MOL/KG)
!       AKVAL(3) - K2C (HCO3 DISSOCIATION) (MOL/KG)
!       AKVAL(4) - KB (B(OH)3 DISSOCIATION) (MOL/KG)  FOR ICONST=4,5,6
!        NKVAL - SIZE OF AKVAL ARRAY (3 FOR ICONST=1,2,3; 4 FOR ICONST=
!        ICALC - SELECTION OF THE TWO INPUT PARAMETERS:
!       ICALC = 1  TOTAL C AND ALKALINITY
!       ICALC = 2  TOTAL C AND PCO2
!       ICALC = 3  TOTAL C AND PH
!       ICALC = 4  ALKALINITY AND PCO2
!       ICALC = 5  ALKALINITY AND PH
!       ICALC = 6  PCO2 AND PH
!       ICALC = 7  CALCULATE CONSTANTS AKVAL ONLY
!       ICONST - SELECTION OF PH SCALE AND COMPONENTS:
!       ICONST = 1  NBS PH SCALE
!       ICONST = 2  HANSSON'S SCALE (SWS WITHOUT FLUORIDE)
!       ICONST = 3  SWS PH SCALE
!       ICONST = 4  AS 1 BUT INCLUDING BORATE IN THE CALCULATION
!       ICONST = 5  AS 2 BUT INCLUDING BORATE IN THE CALCULATION
!       ICONST = 6  AS 3 BUT INCLUDING BORATE IN THE CALCULATION

!  NOTE: FOR ICONST=1,2,3 CONCS(2) REPRESENTS CARBONATE ALKALINITY SINC
!        BORATE IS NOT INCLUDED IN THE CALCULATION. FOR ICONST=4,5,6 CO
!        REPRESENTS TOTAL ALKALINITY (CARBONATE + BORATE), THE COMPONEN
!        WHICH ARE GIVEN IN CONCS(8) AND CONCS(9)
      use fabm_types,only:rk
      use fabm_driver

      IMPLICIT NONE

      real(rk) PMIN, PMAX, SMIN, SMAX, TMIN, TMAX, CONCS,  &
     &      AKVAL, PD, TD, SD, P, T, S, BTOT
      INTEGER MINJC, MAXJC, MINJK, MAXJK, MINCAL, MAXCAL, MINCON,  &
     &      MAXCON, NCONC, NKVAL, ICALC, ICONST, IC
      LOGICAL BORON

      PARAMETER(MINJC=7,MAXJC=9,MINJK=3,MAXJK=4)
      PARAMETER(MINCAL=1,MAXCAL=7,MINCON=1,MAXCON=6)
      PARAMETER(PMIN=0.99999_rk,PMAX=1.00001_rk,SMIN=0.0_rk,  &
     &  SMAX=45.0_rk,TMIN=-4.0_rk,TMAX=40.0_rk)
      DIMENSION CONCS(NCONC),AKVAL(NKVAL)

      P = PD
      S = SD
      T = TD

!     IF(T.LT.TMIN.OR.T.GT.TMAX)WRITE (*,*) P, S, T, TMIN, TMAX
      IF(P.LT.PMIN.OR.P.GT.PMAX) call fatal_error('POLYCO','PRESSURE OUT OF RANGE')
      IF(S.LT.SMIN.OR.S.GT.SMAX) call fatal_error('POLYCO','SALINITY OUT OF RANGE')
      IF(T.LT.TMIN.OR.T.GT.TMAX) call fatal_error('POLYCO','TEMP. OUT OF RANGE')
      IF(ICALC.LT.MINCAL.OR.ICALC.GT.MAXCAL)  &
     &  call fatal_error('POLYCO','ICALC OUT OR RANGE')
      IF(ICONST.LT.MINCON.OR.ICONST.GT.MAXCON)  &
     &  call fatal_error('POLYCO','ICONST OUT OF RANGE')
      BORON=(ICONST.GT.3)
      IF(BORON) THEN
        IC=ICONST-3
        BTOT=0.0004128_rk*S/35.0_rk
        IF(NCONC.NE.MAXJC) call fatal_error('POLYCO','WRONG NCONC VALUE')
        IF(NKVAL.NE.MAXJK) call fatal_error('POLYCO','WRONG NKVAL VALUE')
      ELSE
        IC=ICONST
        IF(NCONC.NE.MINJC) call fatal_error('POLYCO','WRONG NCONC VALUE')
        IF(NKVAL.NE.MINJK) call fatal_error('POLYCO','WRONG NKVAL VALUE')
      ENDIF

      CALL CO2SET(P,T,S,AKVAL,NKVAL,IC)
      IF(ICALC.LT.MAXCAL)  &
     & CALL CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT)

      if (concs(4).eq.100.) then
         write (*,*) 'S,T,P',S,T,P
         write (*,*) 'CONCS',CONCS
         call fatal_error('co2_dyn:POLYCO','Haltafall iteration did not converge.')
      end if

      RETURN
      END SUBROUTINE


      SUBROUTINE CO2SET(P,T,S,AKVAL,NKVAL,IC)
!     -----------------

! Routine to calculate CO2 system constants under the conditions set by
! P,S,T     (NOTE: PRESSURE <> 1ATM IS NOT YET CODED)

! I. Calculate constants at P=1 and S=0 using

!      ln K0  =  A + B/TK + C ln TK
!                                   (where TK is in Kelvin)

! II. Calculate constants at P=1 and salinity S using

!    ln K  =  ln K0 + (a0 + a1/TK + a2 ln TK) S**1/2
!                 + (b0 + b1TK + b2TK**2) S

! The sources of the coefficients are as follows:

!  IC=                  1                    2                  3
!               (NBS pH scale)        (SWS pH scale       (SWS pH scale
!                                       with no F)

!  KP            WEISS (1974)           WEISS(1974)          WEISS(1974

!  K1C )      MEHRBACH ACC. TO       HANSSON ACC. TO    HANSSON AND MEH
!  K2C )       MILLERO (1979)         MILLERO (1979)      ACC. TO DICKS
!   KB )                                                 AND MILLERO (1
!                                                         (K1C AND K2C
!                                                           HANSSON ACC
!                                                            MILLERO (1
!                                                              (KB ONLY

! ***
!      IMPLICIT real(rk) (A-H,O-Z)

!     Modified by jcb 17/02/10 to use OCMIP calculations of K1, K2, Kb.
!     Differences are subtle rather than significant
      use fabm_types

      IMPLICIT NONE

      INTEGER MAXK, MAXCON, NKVAL, ICON, IC, IK
! ***
      PARAMETER(MAXK=4,MAXCON=3)
      real(rk),DIMENSION(:) :: A(MAXK),B(MAXK),C(MAXK)
      real(rk),DIMENSION(:) :: A0(MAXK,MAXCON),A1(MAXK,MAXCON),A2(MAXK,MAXCON)
      real(rk),DIMENSION(:) :: B0(MAXK,MAXCON),B1(MAXK,MAXCON),B2(MAXK,MAXCON)
      real(rk),DIMENSION(:) :: AKVAL(NKVAL)
      real(rk)              :: P,T,S,VAL,TK
      real(rk)              :: dlogTK, S2, sqrtS, S15, k1, k2, kb
      DATA A/-167.8108_rk, 290.9097_rk, 207.6548_rk, 148.0248_rk/
      DATA B/9345.17_rk, -14554.21_rk, -11843.79_rk, -8966.9_rk/
      DATA C/23.3585_rk, -45.0575_rk, -33.6485_rk, -24.4344_rk/
      DATA (A0(1,ICON),ICON=1,MAXCON) /3*0.0_rk/
      DATA (A0(2,ICON),ICON=1,MAXCON) /0.0221_rk, 0.5709_rk, -45.8076_rk/
      DATA (A0(3,ICON),ICON=1,MAXCON) /0.9805_rk, 1.4853_rk, -39.5492_rk/
      DATA (A0(4,ICON),ICON=1,MAXCON) /0.0473_rk, 0.5998_rk, 0.5998_rk/
      DATA (A1(1,ICON),ICON=1,MAXCON) /3*0.0_rk/
      DATA (A1(2,ICON),ICON=1,MAXCON) /34.02_rk, -84.25_rk, 1935.07_rk/
      DATA (A1(3,ICON),ICON=1,MAXCON) /-92.65_rk, -192.69_rk, 1590.14_rk/
      DATA (A1(4,ICON),ICON=1,MAXCON) /49.10_rk, -75.25_rk, -75.25_rk/
      DATA (A2(1,ICON),ICON=1,MAXCON) /3*0.0_rk/
      DATA (A2(2,ICON),ICON=1,MAXCON) /2*0.0_rk,6.9513_rk/
      DATA (A2(3,ICON),ICON=1,MAXCON) /2*0.0_rk,6.1523_rk/
      DATA (A2(4,ICON),ICON=1,MAXCON) /3*0.0_rk/
      DATA (B0(1,ICON),ICON=1,MAXCON) /3*0.023517_rk/
      DATA (B0(2,ICON),ICON=1,MAXCON) /0.0_rk,-0.01632_rk,-0.01566_rk/
      DATA (B0(3,ICON),ICON=1,MAXCON) /-0.03294_rk,-0.05058_rk,-0.04997_rk/
      DATA (B0(4,ICON),ICON=1,MAXCON) /0.0_rk, -0.01767_rk, -0.01767_rk/
      DATA (B1(1,ICON),ICON=1,MAXCON) /3*-2.3656E-4_rk/
      DATA (B1(2,ICON),ICON=1,MAXCON) /3*0.0_rk/
      DATA (B1(3,ICON),ICON=1,MAXCON) /3*0.0_rk/
      DATA (B1(4,ICON),ICON=1,MAXCON) /3*0.0_rk/
      DATA (B2(1,ICON),ICON=1,MAXCON) /3*4.7036E-7_rk/
      DATA (B2(2,ICON),ICON=1,MAXCON) /3*0.0_rk/
      DATA (B2(3,ICON),ICON=1,MAXCON) /3*0.0_rk/
      DATA (B2(4,ICON),ICON=1,MAXCON) /3*0.0_rk/

      TK=T+273.15_rk
      DO 100 IK=1,NKVAL
        VAL=A(IK) + B(IK)/TK + C(IK)*LOG(TK)
        VAL=VAL + (A0(IK,IC) + A1(IK,IC)/TK + A2(IK,IC)*LOG(TK))*SQRT(S)
        VAL=VAL + (B0(IK,IC) + B1(IK,IC)*TK + B2(IK,IC)*TK*TK)*S
        AKVAL(IK)=EXP(VAL)
100    CONTINUE

      IF (IC .EQ. 3) THEN
!  Calculation of constants as used in the OCMIP process for ICONST = 3 or 6
!  see http://www.ipsl.jussieu.fr/OCMIP/
!  added jcb 17/02/10

!  Derive simple terms used more than once
   dlogTK = log(TK)
   S2 = S*S
   sqrtS = sqrt(S)
   S15 = S**1.5
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]
! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
   k1=10**(-1*(3670.7/TK - 62.008 + 9.7944*dlogTK - &
     &      0.0118 * S + 0.000116*S2))
   k2=10**(-1*(1394.7/TK + 4.777 - &
     &      0.0184*S + 0.000118*S2))
! kb = [H][BO2]/[HBO2]
! Millero p.669 (1995) using data from Dickson (1990)
   kb=exp((-8966.90 - 2890.53*sqrtS - 77.942*S + &
     &      1.728*S15 - 0.0996*S2)/TK + &
     &      (148.0248 + 137.1942*sqrtS + 1.62142*S) + &
     &      (-24.4344 - 25.085*sqrtS - 0.2474*S) * &
     &      dlogTK + 0.053105*sqrtS*TK)
! replace haltafall calculations with OCMIP calculations
      AKVAL(2) = k1
      AKVAL(3) = k2
      AKVAL(4) = kb
      END IF ! section implimenting OCMIP coefficients

      RETURN
      END SUBROUTINE


      SUBROUTINE CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT)
!     -----------------

! ROUTINE TO CARRY OUT CO2 CALCULATIONS WITH 2 FIXED PARAMETERS ACCORDI
! THE EQUATIONS GIVEN BY PARKS(1969) AND SKIRROW (1975)
! WITH ADDITIONS FOR INCLUDING BORON IF BORON=.TRUE.


!      IMPLICIT real(rk) (A-H,O-Z)
      use fabm_types

      IMPLICIT NONE

      INTEGER NCONC, NKVAL, ICALC, II, KARL, LQ
      real(rk)              :: CTOT,ALK,PCO2,PH,H2CO3,HCO3,CO3,ALKC
      real(rk)              :: ALKB,AKP,AK1C,AK2C,AKB,BTOT
      real(rk)              :: AKR,AHPLUS
      real(rk)              :: PROD,tol1,tol2,tol3,tol4,steg,fak
      real(rk)              :: STEGBY,Y,X,W,X1,Y1,X2,Y2,FACTOR,TERM,Z
      real(rk),DIMENSION(:) :: CONCS(NCONC),AKVAL(NKVAL),CONCS2(9),AKVAL2(4)
      EQUIVALENCE (CTOT  , CONCS2(1)), (ALK   , CONCS2(2)),  &
     &            (PCO2  , CONCS2(3)), (PH    , CONCS2(4)),  &
     &            (H2CO3 , CONCS2(5)), (HCO3  , CONCS2(6)),  &
     &            (CO3   , CONCS2(7)), (ALKC  , CONCS2(8)),  &
     &            (ALKB  , CONCS2(9)),                       &
     &            (AKP   , AKVAL2(1)), (AK1C  , AKVAL2(2)),  &
     &            (AK2C  , AKVAL2(3)), (AKB   , AKVAL2(4))
      LOGICAL BORON,DONE
      integer :: it
      integer,parameter :: maxiter = 1000


      DO 100 II=1,NCONC
        CONCS2(II)=CONCS(II)
100    CONTINUE
      DO 110 II=1,NKVAL
        AKVAL2(II)=AKVAL(II)
110    CONTINUE
      AKR = AK1C/AK2C
      AHPLUS=10.0_rk**(-PH)
      PROD=AKR*AKP*PCO2

      IF(BORON) THEN

        IF(ICALC.EQ.1.OR.ICALC.EQ.4) THEN
!         *** ALK, BTOT AND CTOT OR PCO2 FIXED ***
!         *** ITERATIVE CALCULATION NECESSARY HERE

!         SET INITIAL GUESSES AND TOLERANCE
          H2CO3=PCO2*AKP
          CO3=ALK/10.0_rk
          AHPLUS=1.0D-8
          ALKB=BTOT
          TOL1=ALK/1.0D5
          TOL2=H2CO3/1.0D5
          TOL3=CTOT/1.0D5
          TOL4=BTOT/1.0D5

!         HALTAFALL iteration to determine CO3, ALKB, AHPLUS
          KARL=1
          STEG=2.0_rk
          FAK=1.0_rk
          STEGBY=0.4_rk
          do it=1,maxiter
             DONE=.TRUE.
             IF(ICALC.EQ.4) THEN
   !         *** PCO2 IS FIXED ***
               Y=AHPLUS*AHPLUS*CO3/(AK1C*AK2C)
               IF(ABS(Y-H2CO3).GT.TOL2) THEN
                 CO3=CO3*H2CO3/Y
                 DONE=.FALSE.
               ENDIF
             ELSEIF(ICALC.EQ.1) THEN
   !           *** CTOT IS FIXED ***
               Y=CO3*(1.0_rk+AHPLUS/AK2C+AHPLUS*AHPLUS/(AK1C*AK2C))
               IF(ABS(Y-CTOT).GT.TOL3) THEN
                 CO3=CO3*CTOT/Y
                 DONE=.FALSE.
               ENDIF
             ENDIF
             Y=ALKB*(1.0_rk+AHPLUS/AKB)
             IF(ABS(Y-BTOT).GT.TOL4) THEN
               ALKB=ALKB*BTOT/Y
               DONE=.FALSE.
             ENDIF

   ! Alkalinity is equivalent to -(total H+), so the sign of W is opposite
   ! to that normally used

             Y=CO3*(2.0_rk+AHPLUS/AK2C)+ALKB
             IF(ABS(Y-ALK).GT.TOL1) THEN
               DONE=.FALSE.
               X=LOG(AHPLUS)
               W=SIGN(1.0_rk,Y-ALK)
               IF(W.GE.0.0_rk) THEN
                 X1=X
                 Y1=Y
               ELSE
                 X2=X
                 Y2=Y
               ENDIF
               LQ=KARL
               IF(LQ.EQ.1) THEN
                 KARL=2*NINT(W)
               ELSEIF(IABS(LQ).EQ.2.AND.(LQ*W).LT.0.) THEN
                 FAK=0.5_rk
                 KARL=3
               ENDIF
               IF(KARL.EQ.3.AND.STEG.LT.STEGBY) THEN
                 W=(X2-X1)/(Y2-Y1)
                 X=X1+W*(ALK-Y1)
               ELSE
                 STEG=STEG*FAK
                 X=X+STEG*W
               ENDIF
               AHPLUS=EXP(X)
             ENDIF
             IF(DONE) exit
          end do
          if (.not. DONE) then
            PH = 100.
          else
             HCO3=CO3*AHPLUS/AK2C
             IF(ICALC.EQ.4) THEN
               CTOT=H2CO3+HCO3+CO3
             ELSEIF(ICALC.EQ.1) THEN
               H2CO3=HCO3*AHPLUS/AK1C
               PCO2=H2CO3/AKP
             ENDIF
             PH=-LOG10(AHPLUS)
             ALKC=ALK-ALKB
          end if
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT, PCO2, AND BTOT FIXED ***
          Y=SQRT(PROD*(PROD-4.0_rk*AKP*PCO2+4.0_rk*CTOT))
          H2CO3=PCO2*AKP
          HCO3=(Y-PROD)/2.0_rk
          CO3=CTOT-H2CO3-HCO3
          ALKC=HCO3+2.0_rk*CO3
          AHPLUS=AK1C*H2CO3/HCO3
          PH=-LOG10(AHPLUS)
          ALKB=BTOT/(1.0_rk+AHPLUS/AKB)
          ALK=ALKC+ALKB
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT, PH AND BTOT FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C)
          CO3=FACTOR*AK1C*AK2C
          HCO3=FACTOR*AK1C*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/AKP
          ALKC=HCO3+2.0_rk*CO3
          ALKB=BTOT/(1.0_rk+AHPLUS/AKB)
          ALK=ALKC+ALKB
        ELSEIF(ICALC.EQ.5) THEN
!         *** ALK, PH AND BTOT FIXED ***
          ALKB=BTOT/(1.0_rk+AHPLUS/AKB)
          ALKC=ALK-ALKB
          HCO3=ALKC/(1.0_rk+2.0_rk*AK2C/AHPLUS)
          CO3=HCO3*AK2C/AHPLUS
          H2CO3=HCO3*AHPLUS/AK1C
          PCO2=H2CO3/AKP
          CTOT=H2CO3+HCO3+CO3
        ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2, PH AND BTOT FIXED ***
          ALKB=BTOT/(1.0_rk+AHPLUS/AKB)
          H2CO3=PCO2*AKP
          HCO3=H2CO3*AK1C/AHPLUS
          CO3=HCO3*AK2C/AHPLUS
          CTOT=H2CO3+HCO3+CO3
          ALKC=HCO3+2.0_rk*CO3
          ALK=ALKC+ALKB
        ENDIF
      ELSE
        IF(ICALC.EQ.1) THEN
!         *** CTOT AND ALK FIXED ***
          TERM=4.0_rk*ALK+CTOT*AKR-ALK*AKR
          Z=SQRT(TERM*TERM+4.0_rk*(AKR-4.0_rk)*ALK*ALK)
          CO3=(ALK*AKR-CTOT*AKR-4.0_rk*ALK+Z)/(2.0_rk*(AKR-4.0_rk))
          HCO3=(CTOT*AKR-Z)/(AKR-4.0_rk)
          H2CO3=CTOT-ALK+CO3
          PCO2=H2CO3/AKP
          PH=-LOG10(AK1C*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT AND PCO2 FIXED ***
          Y=SQRT(PROD*(PROD-4.0_rk*AKP*PCO2+4.0_rk*CTOT))
          H2CO3=PCO2*AKP
          HCO3=(Y-PROD)/2.0_rk
          CO3=CTOT-H2CO3-HCO3
          ALK=HCO3+2.0_rk*CO3
          PH=-LOG10(AK1C*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT AND PH FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C)
          CO3=FACTOR*AK1C*AK2C
          HCO3=FACTOR*AK1C*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/AKP
          ALK=HCO3+2.0_rk*CO3
        ELSEIF(ICALC.EQ.4) THEN
!         *** ALK AND PCO2 FIXED ***
          TERM=SQRT((8.0_rk*ALK+PROD)*PROD)
          CO3=ALK/2.0_rk+PROD/8.0_rk-TERM/8.0_rk
          HCO3=-PROD/4.0_rk+TERM/4.0_rk
          H2CO3=PCO2*AKP
          CTOT=CO3+HCO3+H2CO3
          PH=-LOG10(AK1C*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.5) THEN
!         *** ALK AND PH FIXED ***
          HCO3=ALK/(1.0_rk+2.0_rk*AK2C/AHPLUS)
          CO3=HCO3*AK2C/AHPLUS
          H2CO3=HCO3*AHPLUS/AK1C
          PCO2=H2CO3/AKP
          CTOT=H2CO3+HCO3+CO3
        ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2 AND PH FIXED ***
          H2CO3=PCO2*AKP
          HCO3=H2CO3*AK1C/AHPLUS
          CO3=HCO3*AK2C/AHPLUS
          CTOT=H2CO3+HCO3+CO3
          ALK=HCO3+2.0_rk*CO3
        ENDIF
      ENDIF


      DO 120 II=1,NCONC
        CONCS(II)=CONCS2(II)
120    CONTINUE
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
      SUBROUTINE CaCO3_Saturation (Tc, S, D, CO3, Om_cal, Om_arg)

! Routine to calculate the saturation state of calcite and aragonite
! Inputs:
!               Tc      Temperature (C)
!               S       Salinity
!               D       Depth (m)
!               CO3     Carbonate ion concentration (mol.kg-1 ie /1D6)
!
! Outputs
!               Om_cal  Calite saturation
!               Om_arg  Aragonite saturation
!
! Intermediates
!               K_cal   Stoichiometric solubility product for calcite
!               K_arg   Stoichiometric solubility product for aragonite
!               Ca      Calcium 2+ concentration (mol.kg-1)
!               P       Pressure (bars)
!
! Source
!       Zeebe & Wolf-Gladrow 2001 following Mucci (1983)
!       with pressure corrections from Millero (1995)
!       Code tested against reference values given in Z & W-G
!       Built Jerry Blackford, 2008
!
        use fabm_types

        IMPLICIT None
        real(rk) Tc, Tk, Kelvin, S, D, Ca, CO3
        real(rk) logKspc, Kspc, Om_cal
        real(rk) logKspa, Kspa, Om_arg
        real(rk) tmp1, tmp2, tmp3
        real(rk) dV, dK, P, R

! setup
        Kelvin = 273.15
        Tk = Tc + Kelvin
        Ca = 0.01028    ! Currently oceanic mean value at S=25, needs refining)
        R = 83.131      !(cm3.bar.mol-1.K-1)
        P = D / 10.0    !pressure in bars

! calculate K for calcite
        tmp1 = -171.9065 - (0.077993*Tk) + (2839.319/Tk) + 71.595*log10(Tk)
        tmp2 = + (-0.77712 + (0.0028426*Tk) + (178.34/Tk))*SQRT(S)
        tmp3 = - (0.07711*S) + (0.0041249*(S**1.5))
        logKspc = tmp1 + tmp2 + tmp3
        Kspc = 10.0**logKspc

! correction for pressure for calcite
        IF ( D .GT. 0) THEN
          dV = -48.76 + 0.5304*Tc
          dK = -11.76/1.0D3 + (0.3692/1.0D3) * Tc
          tmp1 = -(dV/(R*Tk))*P + (0.5*dK/(R*Tk))*P*P
          Kspc = Kspc*exp(tmp1)
          logKspc = log10(Kspc)
        END IF

! calculate K for aragonite
        tmp1 = -171.945 - 0.077993*Tk + 2903.293 / Tk + 71.595* log10(Tk)
        tmp2 = + (-0.068393 + 0.0017276*Tk + 88.135/Tk)*SQRT(S)
        tmp3 = - 0.10018*S + 0.0059415*S**1.5
        logKspa = tmp1 + tmp2 + tmp3
        Kspa = 10.0**logKspa

! correction for pressure for aragonite
        IF ( D .GT. 0) THEN
          dV = -46.00 + 0.5304*Tc
          dK = -11.76/1.0D3 + (0.3692/1.0D3) * Tc
          tmp1 = -(dV/(R*Tk))*P + (0.5*dK/(R*Tk))*P*P
          Kspa = Kspa*exp(tmp1)
          logKspa = log10(Kspa)
        END IF

! calculate saturation states
        Om_cal = (CO3 * Ca) / Kspc
        Om_arg = (CO3 * Ca) / Kspa

      RETURN
      END SUBROUTINE

! This file contains a set of FORTRAN subroutines that calculate the carbonate system 
! at any given point in marine space time, given values for 
! temperature, salinity, DIC, depth (pressure). 
! This is essentially an implimentation of the Haltafall speciation code
! (Ingri et al 1967, Talanta 14, 1261 - if it ain't broke don't fix it)
! Another routine calulates the air sea exchange of CO2 given wind speed and atmospheric pCO2.
! Code developed by Jerry blackford and others at PML, based on pre-existing code.
! We accept no liability for errors or inaccuracies, (but will take credit for accuracy and utility!).
! See Zeebe & Wolf-Gladrow, 2001. CO2 in seawater: equilibrium, kinetics and isotopes. 
! Elsevier Oceanography Series 65, 346. for a reasonable overview. 
! Many other packages exist, replicating the same functionality in different languages.
! See http://cdiac.ornl.gov/oceans/co2rprt.html (CO2sys)
! or  http://neon.otago.ac.nz/research/mfc/people/keith_hunter/software/swco2/

! reference for prior usage of this code: Blackford & Gilbert, 2007. J Mar Sys 64, 229-241.

! The following program provides an example implimentation of the routines.
!   PROGRAM Run_CO2_Dynamics
!
!   Real*8 :: Temp, Sal, Depth, DIC, pco2w, pco2a, TA, ph
!   Real*8 :: cco2, carba,bicarb,carb,henry,om_cal,om_arg
!   Real*8 :: wnd, flux
!
!! set up inputs
!   Temp  = 10.0
!   Sal   = 35.0
!   Depth = 0.0
!   DIC   = 2100.0
!   Wnd   = 10.0
!   pCO2a = 390.0
!!!!! note you need to define total alkalinity in subroutine co2_dynamics (see lines 71-86) !!!
!
!   Call CO2_dynamics (Temp, Sal, Depth, DIC,pco2w,TA,ph,cco2, carba,bicarb,carb,henry,om_cal,om_arg)
!
!   Call Air_sea_exchange (Temp, Wnd, pCO2w, pCO2a, Henry, flux)
!
!! write inputs
!   WRITE(*,*) " "
!   WRITE(*,'(A28)') "    .........Inputs........."
!   WRITE(*,'(A18,F10.3)') "    Temperature = ",Temp
!   WRITE(*,'(A18,F10.3)') "    Salinity    = ",Sal
!   WRITE(*,'(A18,F10.3)') "    Depth       = ",Depth
!   WRITE(*,'(A18,F10.3)') "    DIC         = ",DIC
!   WRITE(*,'(A18,F10.3)') "    Wind speed  = ",Wnd
!   WRITE(*,'(A18,F10.3)') "    pCO2 atmos  = ",pCO2a
! 
!! write outputs
!   WRITE(*,*) " "
!   WRITE(*,'(A32,F10.3)') "    ..........Outputs..........."
!   WRITE(*,'(A22,F10.3)') "    pH              = ",pH
!   WRITE(*,'(A22,F10.3)') "    TA              = ",TA
!   WRITE(*,'(A22,F10.3)') "    pco2w           = ",pco2w
!   WRITE(*,'(A22,F10.3)') "    carbonic acid   = ",carba
!   WRITE(*,'(A22,F10.3)') "    bicarbonate     = ",bicarb
!   WRITE(*,'(A22,F10.3)') "    carbonate       = ",carb
!   WRITE(*,'(A22,F10.3)') "    Omega calcite   = ",om_cal
!   WRITE(*,'(A22,F10.3)') "    Omega aragonite = ",om_arg
!   WRITE(*,'(A22,F10.3)') "    air sea flux    = ",flux
!   WRITE(*,*) " "
!
!
!   STOP
!   END PROGRAM Run_CO2_Dynamics




   subroutine CO2_dynamics(T,S,Z,DIC,pco2w,TA,ph,cco2, carba,bicarb,carb,henry,om_cal,om_arg)
  
   IMPLICIT NONE

! INPUT PARAMETERS:
      real*8   :: T, S, Z, DIC
! T: temperature (C)
! S: salinity (psu)
! Z: depth (metres)
! DIC: total dissolved inorganic carbon (also commonly called TC total carbon)
! pCO2a: partial pressure of co2 in atmosphere (eg ~390 today, 280=preindustrial, sky's the limit)

! LOCAL VARIABLES:
      real*8   :: TCO2, TA, ph, pco2w
      real*8   :: pco2water,fairco2, pco2x, henry
      real*8   :: ca, bc, cb, cco2, carba, bicarb, carb, om_cal, om_arg

!!!!!!!! provide an expression for total alkalinity.!!!!!!!!!!!!!!
! Note, for oceanic regimes there is generally a well constrained relationship between salinity and total alkalinity (tA) as TA is conservative. However the slope and intercept vary according to region.
! For coastal regions, riverine input of TA makes the TA - s relationship unsatisfactory.
! If you have the data then one possibility is to introduce riverine alkalinity as a seperate tracer.
! note that uptake/release of nutrients and Calcium carbonate also make small modifications to local TA.
! See Zeebe & Wolf-Gladrow, 
!     Millero et al 1998, Mar Chem 60, 111-130. 
!     Lee et al, 2006, GRL 33, L19605, 
!     Wolf-Gladrow et al, Mar Chem 106 (2007).

! some example relationships with salinity

!          TA = 66.96*S - 36.803 ! Nordic, N Atlantic, Bellerby et al 2005. 
          TA = 520.1 + 51.24*S  ! Atlantic (Millero 1998)
!          TA = 399.0 + 54.629*S ! Pacific (Millero, 1998)
!          TA = -114.0 + 68.80*S ! Indian (Millero, 1998)

!scale inputs for iteration (this code uses mol/kg, inputs are generally in mmol/m3)
          TA = TA / 1.0D6
          TCO2  = DIC / 1.0D6

!..call the parent routine for the carbonate system
          CALL CO2DYN ( TCO2, TA, T, S, PCO2WATER, pH, HENRY, PCO2X, ca, bc, cb)

!scale the outputs of the routine to useful units (generally mmols C m-3)and name as appropriate 
          PCO2W = PCO2WATER*1.0D6 ! partial pressure of co2 in water
          TA = TA*1.0D6         ! total alkalinity
          Cco2 = PCO2X * 1.0D6    ! concentration of co2 (partial pressure x solubility)
          CarbA= ca*1.0d6         ! carbonic acid concentration
          Bicarb=bc*1.0d6         ! bicarbonate ion concentration
          Carb = cb*1.0d6         ! carbonate ion concentration
!         pH                         ! pH, correct units
!         Henry                      ! Henry's constant, correct units

!Call carbonate saturation state subroutine to calculate calcite and aragonite calcification states
          CALL CaCO3_Saturation ( T, S, Z, cb, Om_cal, Om_arg)
! outputs are
!         Om_cal calcite saturation state
!         Om_arg aragonite saturation state (both ok units)

   RETURN

   END SUBROUTINE CO2_dynamics


   SUBROUTINE Air_sea_exchange (T, Wnd, pCO2w, pCO2a, Henry, flux)

!  this routine should be called for the surface box only 
!  Uses the Nightingale and Liss parameterisation (GBC 14, 373-388 (2000).
!  it calculates the air-sea flux of CO2 given
!  pCO2w    partial pressure of CO2 in the water (from carbonate system subroutine call)
!  pCO2a    partial pressure of CO2 in the atmosphere (usually external forcing).
!  T        temperature (C)
!  Wnd      wind speed, metres
!  Henry    henry's constant
!  Outputs are
!  flux     flux of CO2 in mmol C /m2/d +ve is in-gassing (air to sea), -ve is outgassing (sea to air).
!      SC is the Schmidt number from Wanninkhof (1992, JGR 97,7373-7382), the Nightingale et al (2000)
!      transfer velocity being for Sc=600 and the dependence of tranfer velocity on Sc coming 
!      from Jahne et al (1987, J. Geophys Res 92, 1937-1949).

   IMPLICIT NONE

! INPUT PARAMETERS:
      real*8   :: T, wnd, pco2w, pco2a, henry
! LOCAL VARIABLES:
      real*8   :: sc, fwind
! OUTPUT Variables
      real*8   :: flux

! calculate the scmidt number and unit conversions
          sc=2073.1-125.62*T+3.6276*T**2.0-0.0432190*T**3.0
          fwind =  (0.222d0 * wnd**2d0 + 0.333d0 * wnd)*(sc/660.d0)**(-0.5)
          fwind=fwind*24.d0/100.d0   ! convert to m/day

! flux depends on the difference in partial pressures, wind and henry
          flux = fwind * HENRY * ( PCO2A - PCO2W )

  RETURN 
  END SUBROUTINE Air_sea_exchange



      SUBROUTINE CO2dyn ( TCO2, TA, T, S, pco2, ph, henry, pco2x, ca, bc, cb)

!.......................................................................
!     This subroutine acts as an interface to the Haltafall code, defining parameters etc.
!.......................................................................


      real*8 PRSS, PH, AKVAL, CONCS,                       &
     &                 TCO2, TA, T, S, PCO2, PCO2X,           & 
     &                 SOLBTY, CCO2,                                 &
     &                 A1, A2, A3, B1, B2, B3, TK, TK1, SOL1, SOL2,  &
     &                 HENRY, ca, bc, cb
      INTEGER MCONC, MKVAL, ICONST, ICALC

! INPUT PARAMETERS:
      PARAMETER ( MCONC = 9,MKVAL = 4 )

      DIMENSION AKVAL(MKVAL), CONCS(MCONC)

      ICONST = 6
      PRSS = 1.0D0
      CONCS(1) = TCO2
      CONCS(2) = TA
      ICALC = 1
      
      CALL POLYCO(PRSS,T,S,CONCS,MCONC,AKVAL,MKVAL,ICALC,ICONST)

      PCO2   = CONCS(3)
      PH = CONCS(4)
      ca = CONCS(5)
      bc = CONCS(6)
      cb = CONCS(7)

!..calculate solubility of CO2........................................
      TK     = T + 273.15D0

      A1 = 0.0D0 - 60.2409D0
      A2 = 93.4517D0
      A3 = 23.3585D0
      B1 = 0.023517D0
      B2 = 0.0D0 - 0.023656D0
      B3 = 0.0047036D0
      TK1    = TK / 100.0D0
      SOL1   = A1 + (A2/TK1) + (A3*LOG(TK1))
      SOL2   = 36.4D0 * (B1 + (B2*TK1) + (B3*TK1*TK1))
      SOLBTY = EXP((SOL1 + SOL2))

!     concentration of co2 = partial pressure of co2 * solubility
      PCO2X   = TCO2 / SOLBTY

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
!       ICONST = 5  AS 3 BUT INCLUDING BORATE IN THE CALCULATION

!  NOTE: FOR ICONST=1,2,3 CONCS(2) REPRESENTS CARBONATE ALKALINITY SINC
!        BORATE IS NOT INCLUDED IN THE CALCULATION. FOR ICONST=4,5,6 CO
!        REPRESENTS TOTAL ALKALINITY (CARBONATE + BORATE), THE COMPONEN
!        WHICH ARE GIVEN IN CONCS(8) AND CONCS(9)


      real*8 PMIN, PMAX, SMIN, SMAX, TMIN, TMAX, CONCS,  &
     &      AKVAL, PD, TD, SD, P, T, S, BTOT
      INTEGER MINJC, MAXJC, MINJK, MAXJK, MINCAL, MAXCAL, MINCON,  &
     &      MAXCON, NCONC, NKVAL, ICALC, ICONST, IC
      LOGICAL BORON

      PARAMETER(MINJC=7,MAXJC=9,MINJK=3,MAXJK=4)
      PARAMETER(MINCAL=1,MAXCAL=7,MINCON=1,MAXCON=6)
      PARAMETER(PMIN=0.99999D0,PMAX=1.00001D0,SMIN=0.0D0,  &
     &  SMAX=40.0D0,TMIN=0.0D0,TMAX=40.0D0)
      DIMENSION CONCS(NCONC),AKVAL(NKVAL)

      P = PD
      S = SD
      T = TD

!     IF(T.LT.TMIN.OR.T.GT.TMAX)WRITE (*,*) P, S, T, TMIN, TMAX
      IF(P.LT.PMIN.OR.P.GT.PMAX) STOP'POLYCO - PRESSURE OUT OF RANGE'
      IF(S.LT.SMIN.OR.S.GT.SMAX) STOP'POLYCO - SALINITY OUT OF RANGE'
      IF(T.LT.TMIN.OR.T.GT.TMAX) STOP'POLYCO - TEMP. OUT OF RANGE'
      IF(ICALC.LT.MINCAL.OR.ICALC.GT.MAXCAL)  &
     &  STOP'POLYCO - ICALC OUT OR RANGE'
      IF(ICONST.LT.MINCON.OR.ICONST.GT.MAXCON)  &
     &  STOP'POLYCO - ICONST OUT OF RANGE'
      BORON=(ICONST.GT.3)
      IF(BORON) THEN
        IC=ICONST-3
        BTOT=0.0004128D0*S/35.0D0
        IF(NCONC.NE.MAXJC) STOP'POLYCO - WRONG NCONC VALUE'
        IF(NKVAL.NE.MAXJK) STOP'POLYCO - WRONG NKVAL VALUE'
      ELSE
        IC=ICONST
        IF(NCONC.NE.MINJC) STOP'POLYCO - WRONG NCONC VALUE'
        IF(NKVAL.NE.MINJK) STOP'POLYCO - WRONG NKVAL VALUE'
      ENDIF

      CALL CO2SET(P,T,S,AKVAL,NKVAL,IC)
      IF(ICALC.LT.MAXCAL)  &
     & CALL CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT)
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
!      IMPLICIT real*8 (A-H,O-Z)
      INTEGER MAXK, MAXCON, NKVAL, ICON, IC, IK
! ***
      PARAMETER(MAXK=4,MAXCON=3)
      real*8,DIMENSION(:) :: A(MAXK),B(MAXK),C(MAXK)
      real*8,DIMENSION(:) :: A0(MAXK,MAXCON),A1(MAXK,MAXCON),A2(MAXK,MAXCON)
      real*8,DIMENSION(:) :: B0(MAXK,MAXCON),B1(MAXK,MAXCON),B2(MAXK,MAXCON)
      real*8,DIMENSION(:) :: AKVAL(NKVAL)
      real*8              :: P,T,S,VAL,TK
      DATA A/-167.8108D0, 290.9097D0, 207.6548D0, 148.0248D0/
      DATA B/9345.17D0, -14554.21D0, -11843.79D0, -8966.9D0/
      DATA C/23.3585D0, -45.0575D0, -33.6485D0, -24.4344D0/
      DATA (A0(1,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (A0(2,ICON),ICON=1,MAXCON) /0.0221D0, 0.5709D0, -45.8076D0/
      DATA (A0(3,ICON),ICON=1,MAXCON) /0.9805D0, 1.4853D0, -39.5492D0/
      DATA (A0(4,ICON),ICON=1,MAXCON) /0.0473D0, 0.5998D0, 0.5998D0/
      DATA (A1(1,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (A1(2,ICON),ICON=1,MAXCON) /34.02D0, -84.25D0, 1935.07D0/
      DATA (A1(3,ICON),ICON=1,MAXCON) /-92.65D0, -192.69D0, 1590.14D0/
      DATA (A1(4,ICON),ICON=1,MAXCON) /49.10D0, -75.25D0, -75.25D0/
      DATA (A2(1,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (A2(2,ICON),ICON=1,MAXCON) /2*0.0D0,6.9513D0/
      DATA (A2(3,ICON),ICON=1,MAXCON) /2*0.0D0,6.1523D0/
      DATA (A2(4,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B0(1,ICON),ICON=1,MAXCON) /3*0.023517D0/
      DATA (B0(2,ICON),ICON=1,MAXCON) /0.0D0,-0.01632D0,-0.01566D0/
      DATA (B0(3,ICON),ICON=1,MAXCON) /-0.03294D0,-0.05058D0,-0.04997D0/
      DATA (B0(4,ICON),ICON=1,MAXCON) /0.0D0, -0.01767D0, -0.01767D0/
      DATA (B1(1,ICON),ICON=1,MAXCON) /3*-2.3656D-4/
      DATA (B1(2,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B1(3,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B1(4,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B2(1,ICON),ICON=1,MAXCON) /3*4.7036D-7/
      DATA (B2(2,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B2(3,ICON),ICON=1,MAXCON) /3*0.0D0/
      DATA (B2(4,ICON),ICON=1,MAXCON) /3*0.0D0/

      TK=T+273.15D0
      DO 100 IK=1,NKVAL
        VAL=A(IK) + B(IK)/TK + C(IK)*LOG(TK)
        VAL=VAL + (A0(IK,IC) + A1(IK,IC)/TK + A2(IK,IC)*LOG(TK))*DSQRT(S)
        VAL=VAL + (B0(IK,IC) + B1(IK,IC)*TK + B2(IK,IC)*TK*TK)*S
        AKVAL(IK)=EXP(VAL)
100    CONTINUE

      RETURN
      END SUBROUTINE


      SUBROUTINE CO2CLC(CONCS,NCONC,AKVAL,NKVAL,ICALC,BORON,BTOT)
!     -----------------

! ROUTINE TO CARRY OUT CO2 CALCULATIONS WITH 2 FIXED PARAMETERS ACCORDI
! THE EQUATIONS GIVEN BY PARKS(1969) AND SKIRROW (1975)
! WITH ADDITIONS FOR INCLUDING BORON IF BORON=.TRUE.


!      IMPLICIT real*8 (A-H,O-Z)
      INTEGER NCONC, NKVAL, ICALC, II, KARL, LQ
      real*8              :: CTOT,ALK,PCO2,PH,H2CO3,HCO3,CO3,ALKC
      real*8              :: ALKB,AKP,AK1C,AK2C,AKB,BTOT
      real*8              :: AKR,AHPLUS
      real*8              :: PROD,tol1,tol2,tol3,tol4,steg,fak
      real*8              :: STEGBY,Y,X,W,X1,Y1,X2,Y2,FACTOR,TERM,Z
      real*8,DIMENSION(:) :: CONCS(NCONC),AKVAL(NKVAL),CONCS2(9),AKVAL2(4)
      EQUIVALENCE (CTOT  , CONCS2(1)), (ALK   , CONCS2(2)),  &
     &            (PCO2  , CONCS2(3)), (PH    , CONCS2(4)),  &
     &            (H2CO3 , CONCS2(5)), (HCO3  , CONCS2(6)),  &
     &            (CO3   , CONCS2(7)), (ALKC  , CONCS2(8)),  &
     &            (ALKB  , CONCS2(9)),                       &
     &            (AKP   , AKVAL2(1)), (AK1C  , AKVAL2(2)),  &
     &            (AK2C  , AKVAL2(3)), (AKB   , AKVAL2(4))
      LOGICAL BORON,DONE


      DO 100 II=1,NCONC
        CONCS2(II)=CONCS(II)
100    CONTINUE
      DO 110 II=1,NKVAL
        AKVAL2(II)=AKVAL(II)
110    CONTINUE
      AKR = AK1C/AK2C
      AHPLUS=10.0D0**(-PH)
      PROD=AKR*AKP*PCO2
     
      IF(BORON) THEN

        IF(ICALC.EQ.1.OR.ICALC.EQ.4) THEN
!         *** ALK, BTOT AND CTOT OR PCO2 FIXED ***
!         *** ITERATIVE CALCULATION NECESSARY HERE

!         SET INITIAL GUESSES AND TOLERANCE
          H2CO3=PCO2*AKP
          CO3=ALK/10.0D0
          AHPLUS=1.0D-8
          ALKB=BTOT
          TOL1=ALK/1.0D5
          TOL2=H2CO3/1.0D5
          TOL3=CTOT/1.0D5
          TOL4=BTOT/1.0D5

!         HALTAFALL iteration to determine CO3, ALKB, AHPLUS
          KARL=1
          STEG=2.0D0
          FAK=1.0D0
          STEGBY=0.4D0
10        DONE=.TRUE.
          IF(ICALC.EQ.4) THEN
!         *** PCO2 IS FIXED ***
            Y=AHPLUS*AHPLUS*CO3/(AK1C*AK2C)
            IF(ABS(Y-H2CO3).GT.TOL2) THEN
              CO3=CO3*H2CO3/Y
              DONE=.FALSE.
            ENDIF
          ELSEIF(ICALC.EQ.1) THEN
!           *** CTOT IS FIXED ***
            Y=CO3*(1.0D0+AHPLUS/AK2C+AHPLUS*AHPLUS/(AK1C*AK2C))
            IF(ABS(Y-CTOT).GT.TOL3) THEN
              CO3=CO3*CTOT/Y
              DONE=.FALSE.
            ENDIF
          ENDIF
          Y=ALKB*(1.0D0+AHPLUS/AKB)
          IF(ABS(Y-BTOT).GT.TOL4) THEN
            ALKB=ALKB*BTOT/Y
            DONE=.FALSE.
          ENDIF

! Alkalinity is equivalent to -(total H+), so the sign of W is opposite
! to that normally used

          Y=CO3*(2.0D0+AHPLUS/AK2C)+ALKB
          IF(ABS(Y-ALK).GT.TOL1) THEN
            DONE=.FALSE.
            X=LOG(AHPLUS)
            W=SIGN(1.0D0,Y-ALK)
            IF(W.GE.0.0D0) THEN
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
              FAK=0.5D0
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
          IF(.NOT.DONE) GOTO 10
          
          HCO3=CO3*AHPLUS/AK2C
          IF(ICALC.EQ.4) THEN
            CTOT=H2CO3+HCO3+CO3
          ELSEIF(ICALC.EQ.1) THEN
            H2CO3=HCO3*AHPLUS/AK1C
            PCO2=H2CO3/AKP
          ENDIF
          PH=-LOG10(AHPLUS)
          ALKC=ALK-ALKB
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT, PCO2, AND BTOT FIXED ***
          Y=SQRT(PROD*(PROD-4.0D0*AKP*PCO2+4.0D0*CTOT))
          H2CO3=PCO2*AKP
          HCO3=(Y-PROD)/2.0D0
          CO3=CTOT-H2CO3-HCO3
          ALKC=HCO3+2.0D0*CO3
          AHPLUS=AK1C*H2CO3/HCO3
          PH=-LOG10(AHPLUS)
          ALKB=BTOT/(1.0D0+AHPLUS/AKB)
          ALK=ALKC+ALKB
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT, PH AND BTOT FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C)
          CO3=FACTOR*AK1C*AK2C
          HCO3=FACTOR*AK1C*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/AKP
          ALKC=HCO3+2.0D0*CO3
          ALKB=BTOT/(1.0D0+AHPLUS/AKB)
          ALK=ALKC+ALKB
        ELSEIF(ICALC.EQ.5) THEN
!         *** ALK, PH AND BTOT FIXED ***
          ALKB=BTOT/(1.0D0+AHPLUS/AKB)
          ALKC=ALK-ALKB
          HCO3=ALKC/(1.0D0+2.0D0*AK2C/AHPLUS)
          CO3=HCO3*AK2C/AHPLUS
          H2CO3=HCO3*AHPLUS/AK1C
          PCO2=H2CO3/AKP
          CTOT=H2CO3+HCO3+CO3
        ELSEIF(ICALC.EQ.6) THEN
!         *** PCO2, PH AND BTOT FIXED ***
          ALKB=BTOT/(1.0D0+AHPLUS/AKB)
          H2CO3=PCO2*AKP
          HCO3=H2CO3*AK1C/AHPLUS
          CO3=HCO3*AK2C/AHPLUS
          CTOT=H2CO3+HCO3+CO3
          ALKC=HCO3+2.0D0*CO3
          ALK=ALKC+ALKB
        ENDIF
      ELSE
        IF(ICALC.EQ.1) THEN
!         *** CTOT AND ALK FIXED ***
          TERM=4.0D0*ALK+CTOT*AKR-ALK*AKR
          Z=SQRT(TERM*TERM+4.0D0*(AKR-4.0D0)*ALK*ALK)
          CO3=(ALK*AKR-CTOT*AKR-4.0D0*ALK+Z)/(2.0D0*(AKR-4.0D0))
          HCO3=(CTOT*AKR-Z)/(AKR-4.0D0)
          H2CO3=CTOT-ALK+CO3
          PCO2=H2CO3/AKP
          PH=-LOG10(AK1C*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.2) THEN
!         *** CTOT AND PCO2 FIXED ***
          Y=SQRT(PROD*(PROD-4.0D0*AKP*PCO2+4.0D0*CTOT))
          H2CO3=PCO2*AKP
          HCO3=(Y-PROD)/2.0D0
          CO3=CTOT-H2CO3-HCO3
          ALK=HCO3+2.0D0*CO3
          PH=-LOG10(AK1C*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.3) THEN
!         *** CTOT AND PH FIXED ***
          FACTOR=CTOT/(AHPLUS*AHPLUS+AK1C*AHPLUS+AK1C*AK2C)
          CO3=FACTOR*AK1C*AK2C
          HCO3=FACTOR*AK1C*AHPLUS
          H2CO3=FACTOR*AHPLUS*AHPLUS
          PCO2=H2CO3/AKP
          ALK=HCO3+2.0D0*CO3
        ELSEIF(ICALC.EQ.4) THEN
!         *** ALK AND PCO2 FIXED ***
          TERM=SQRT((8.0D0*ALK+PROD)*PROD)
          CO3=ALK/2.0D0+PROD/8.0D0-TERM/8.0D0
          HCO3=-PROD/4.0D0+TERM/4.0D0
          H2CO3=PCO2*AKP
          CTOT=CO3+HCO3+H2CO3
          PH=-LOG10(AK1C*H2CO3/HCO3)
        ELSEIF(ICALC.EQ.5) THEN
!         *** ALK AND PH FIXED ***
          HCO3=ALK/(1.0D0+2.0D0*AK2C/AHPLUS)
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
          ALK=HCO3+2.0D0*CO3
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
!               CO3     Carbonate ion concentration (mol.-l ie /1D6)    
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

        IMPLICIT None
        REAL*8 Tc, Tk, Kelvin, S, D, Ca, CO3
        REAL*8 logKspc, Kspc, Om_cal
        REAL*8 logKspa, Kspa, Om_arg
        REAL*8 tmp1, tmp2, tmp3
        REAL*8 dV, dK, P, R
               
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
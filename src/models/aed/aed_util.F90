!###############################################################################
!#                                                                             #
!# aed_util.F90                                                                #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created June  2012                                                          #
!#                                                                             #
!###############################################################################

#ifdef _FABM_F2003_

#include "aed.h"

!
MODULE aed_util
!-------------------------------------------------------------------------------
! aed_util --- shared utility functions for aed modules
!
!-------------------------------------------------------------------------------
   USE fabm_types

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC find_free_lun, qsort
   PUBLIC aed_gas_piston_velocity, aed_oxygen_sat, exp_integral
   PUBLIC aed_bio_temp_function,fTemp_function
!


!===============================================================================
CONTAINS



!###############################################################################
INTEGER FUNCTION find_free_lun()
!-------------------------------------------------------------------------------
! find a free logical unit number
!-------------------------------------------------------------------------------
!LOCALS
    INTEGER :: lun
    LOGICAL :: opend
!
!-------------------------------------------------------------------------------
!BEGIN
   DO lun = 10,99
      inquire( unit = lun, opened = opend )
      IF ( .not. opend ) THEN
         find_free_lun = lun
         RETURN
      ENDIF
   ENDDO

   find_free_lun = -1
END FUNCTION find_free_lun
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE real(rk) FUNCTION aed_gas_piston_velocity(wshgt,wind,tem,sal)
!-------------------------------------------------------------------------------
! Atmospheric-surface water exchange piston velocity for O2, CO2 etc
!-------------------------------------------------------------------------------
!ARGUMENTS
   real(rk),INTENT(IN)    :: wshgt,wind,tem,sal
!
!LOCALS
   ! Temporary variables
   real(rk) :: schmidt,k_wind,k_flow,temp,salt,hgtCorrx
   ! Parameters
   real(rk),PARAMETER :: roughlength = 0.000114  ! momn roughness length(m)
   real(rk),PARAMETER :: SCHMIDT_MODEL = 2.0
!
!-------------------------------------------------------------------------------
!BEGIN
   temp=tem
   salt=sal
   IF (temp < 0.0)       temp = 0.0
   IF (temp > 38.0)      temp = 38.0
   IF (salt < 0.0)       salt = 0.0
   IF (salt > 75.0)      salt = 75.0

   k_flow = 0.0_rk ! Needs to be set based on flow velocity

   ! Schmidt, Sc
   ! control value : Sc = 590 at 20°C and 35 psu

    schmidt = 590.

    IF(SCHMIDT_MODEL==1) THEN
      schmidt = (0.9 + 0.1*salt/35.0)*(1953.4+temp*(-128.0+temp*(3.9918-temp*0.050091)))
    ELSEIF(SCHMIDT_MODEL==2) THEN
      schmidt = (0.9 + salt/350.0)
      schmidt = schmidt * (2073.1 -125.62*temp +3.6276*temp*temp - 0.043219*temp*temp*temp)
    ELSEIF(SCHMIDT_MODEL==3) THEN
      ! http://www.geo.uu.nl/Research/Geochemistry/kb/Knowledgebook/O2_transfer.pdf
      schmidt = (1.0 + 3.4e-3*salt)
      schmidt = schmidt * (1800.6 -120.1*temp +3.7818*temp*temp - 0.047608*temp*temp*temp)
    ENDIF


   ! Adjust the windspeed if the sensor height is not 10m

    hgtCorrx =  LOG(10.00 / roughLength) / LOG(wshgt / roughLength)

   ! Gas transfer velocity, kCO2 (cm/hr)
   ! k = 0.31 u^2 (Sc/660)^-0.5
   ! This parameterization of course assumes 10m windspeed, and so
   ! must be scaled by hgtCorrx

    k_wind = 0.31 * wind*wind*hgtCorrx*hgtCorrx / SQRT(schmidt/660.0) !in cm/hr

   ! convert to m/s

    k_wind = k_wind / 3.6e5

   ! piston velocity is the sum due to flow and wind

    aed_gas_piston_velocity = k_flow + k_wind
END FUNCTION aed_gas_piston_velocity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE real(rk) FUNCTION aed_oxygen_sat(salt,temp)
!-------------------------------------------------------------------------------
!  Calculated saturated oxygen concentration at salinity and temperature
! Taken from Riley and Skirrow (1974)
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   real(rk),INTENT(in) :: salt,temp
!
!LOCALS
   real(rk) :: Tabs
   real(rk) :: buf1, buf2, buf3, sol_coeff
!
!-------------------------------------------------------------------------------
!BEGIN
   buf1 = 0.0_rk ; buf2 = 0.0_rk ; buf3 = 0.0_rk ; sol_coeff = 0.0_rk

   Tabs = temp + 273.15
   buf1 = -173.4292 + 249.6339 * 100.0 / Tabs + 143.3483 * LOG(Tabs/100.0)
   buf2 = 21.8492 * Tabs / 100.0
   buf3 = salt * (-0.033096 + 0.014259 * Tabs / 100.0 - 0.0017 * (Tabs / 100.0)**2.0)
   sol_coeff = buf1 - buf2 + buf3

   aed_oxygen_sat = 1.42763 * exp(sol_coeff) !in g/m3

   !Convert to mmol/m3
   aed_oxygen_sat = (aed_oxygen_sat / 32.) * 1e3
END FUNCTION aed_oxygen_sat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION exp_integral(inp) RESULT(E_ib)
!ARGUMENTS
   real(rk)  :: inp
!
!LOCALS
   real(rk)  :: E_ib !-- Outgoing
   INTEGER   :: j
   real(rk)  :: ff
!
!-------------------------------------------------------------------------------
!BEGIN
   ff = -1e-9
   IF(ABS(inp-10.0) < 12.0) THEN
     IF(inp==0.0) THEN
       E_ib = inp
     ELSE
       j  = 10+2*IABS(INT(inp))
       ff = 1.0/(REAL(j+1)**2.0)
       DO WHILE(j/=0)
         ff = (ff*REAL(j)*inp+1.0)/REAL(j*j)
         j  = j-1
       ENDDO
       ff   = ff*inp+LOG(1.781072418*ABS(inp))
       E_ib = ff
     ENDIF
   ELSE
     j = 5 + 20 / IABS(INT(inp))
     ff = inp
     DO WHILE(j/=0)
       ff = (1.0/(1.0/ff-1.0/REAL(j)))+inp
       j = j-1
     ENDDO
     ff  = EXP(inp)/ff
     E_ib = ff
   ENDIF

END FUNCTION exp_integral
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_bio_temp_function(numg, theta, T_std, T_opt, T_max, aTn, bTn, kTn, name)
!-------------------------------------------------------------------------------
! Numerically solver for continuos temperature function
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)       :: numg        ! Number of groups
   real(rk),INTENT(in)      :: theta(:)
   real(rk),INTENT(inout)   :: T_std(:), T_opt(:), T_max(:)
   real(rk),INTENT(out)     :: aTn(:), bTn(:), kTn(:)
   CHARACTER(64),INTENT(in) :: name(:)
!
!LOCALS
   real(rk) :: Ts     ! Min. temperature where fT(Ts)=I (usually 1)
   real(rk) :: To     ! Optimum temperature where d(fT(To))/dT=0
   real(rk) :: Tm     ! Maximum temperature where fT(Tm)=0
   real(rk) :: in     ! Constant for fT(Ts)=in
   real(rk) :: v      ! Constant v
   real(rk) :: k,a,b  ! Model constants
   real(rk) :: G      ! Function fT()
   real(rk) :: devG       ! Derivative of fT()
   real(rk) :: a0,a1,a2   ! Dummies
   real(rk) :: tol        ! Tolerance
   INTEGER :: group       ! Group counter
   INTEGER :: i           ! Counters
   real(rk),PARAMETER :: t20=20.0
   LOGICAL,PARAMETER :: curvef=.true. ! T : f(T)=v**(T-20) at T=Tsta
                                      ! F : f(T) = 1 at T=Tsta

   real(rk),ALLOCATABLE,DIMENSION(:,:) :: value
!
!-------------------------------------------------------------------------------
!BEGIN
    write(*,"('Estimating temperature functions for phytoplankton - ')")
    write(*,"(' Temperature function of the form :',/)")
    write(*,"('    fT = v^(T-20)-v^(k(T-a))+b',/)")

    tol   = 0.05

    DO group=1,numg

      ! Set the constants for the correct group
      v = theta(group)

      IF(v < 1.01) THEN
        print "(/,2X,'WARNING: theta_growth for group ',I2,' < 1.01',/)",group
      ENDIF

      Tm = T_max(group)
      Ts = T_std(group)
      To = T_opt(group)


      IF (Ts<0.0 .AND. To<0.0 .AND. Tm<0.0) THEN
        ! The user inputs the values of kTn, aTn and bTn directly
        kTn(group) = -Ts
        bTn(group) = -Tm
        aTn(group) = -To

        ALLOCATE(value(401,1))
        ! Calculate the temperature function using 0.1 deg C intervals

        DO i = 0,400
          b = REAL(i)/10.0
          value(i+1,1) = v**(b-20) - v**(kTn(group) * (b - aTn(group))) + bTn(group)
        ENDDO

        ! Find the values of Tsta, T_opt and T_max from the temp function
        a=0.0
        DO i=1,SIZE(value,1)
          b=REAL(i-1)/10.0
          IF(value(i,1)>0.0) THEN
            T_max(group) = b
          ENDIF
          IF(value(i,1)>a) THEN
            T_opt(group) = b
            a=value(i,1)
          ENDIF
          IF(value(i,1)>v**(b-20)-tol .and. value(i,1)<v**(b-20)+tol) THEN
            T_std(group) = b
          ENDIF
        ENDDO
        DEALLOCATE(value)


      ELSE
        in = 1.0
        a0 = v**(Ts-t20)
        a1 = v**(To-t20)
        a2 = v**(Tm-t20)

        ! Perform the iteration to find the constants.
        ! First approximation of k.
        k = 6.0
        i = 0
        G = tol + 1.0
        ! Do the iterations until -tol < G < tol
        DO WHILE((G <= -tol) .OR. (G >= tol))
          i=i+1
          IF(i==100) THEN  ! Increases the tolerance if more than 100
            i=0            ! iterations are performed.
            tol=tol+0.01
          ENDIF
          IF(curvef) THEN
            ! Use the condition f(T)=v**(T-20) at T=Tsta
            G = k * v**(k * To) * a2 - a1 * (v**(k * Tm) - v**(k * Ts))
            devG = v**(k * To) * a2 * (in + k * To * log(v)) - a1 * log(v) &
              * (Tm * v**(k * Tm) - Ts * v**(k * Ts))
          ELSE
            ! Use the condition f(T)=1 at T=Tsta
            G = k * v**(k * To) * (a0 - a2 - in) - a1 * (v**(k * Ts) &
              - v**(k * Tm))
            devG = (a0 - a2 - in) * v**(k * To) * (in + k * To * log(v)) - a1 &
              * log(v) * (Ts * v**(k * Ts) - Tm * v**(k * Tm))
          ENDIF
          ! Find the next iteration of k
          k = k - G / devG
        ENDDO

        ! Get the remaining model constants
        IF(k/=0.0) THEN
          a=-log(a1/(k*v**(k*To)))/(k*log(v))
          IF(curvef) THEN
            b=v**(k*(Ts-a))
          ELSE
            b=in+v**(k*(Ts-a))-a0
          ENDIF
        ELSE
          a=0.0
          b=0.0
        ENDIF

        ! Set the model constants to the calculated values
        kTn(group) = k
        aTn(group) = a
        bTn(group) = b
      ENDIF

      IF (kTn(group) < 0.1 .AND. bTn(group) > 100.0) THEN
            PRINT *,'Cannot solve for fT for: ', name(group)
            STOP
      ENDIF

    ENDDO

END SUBROUTINE aed_bio_temp_function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION fTemp_function(method,T_max,T_std,theta,aTn,bTn,kTn,temp) RESULT(fT)
!-------------------------------------------------------------------------------
! Generic temperature function for phytoplankton and zooplankton taking into
! account a decrease in production above T_opt.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)  :: method
   real(rk),INTENT(in) :: T_max, T_std,theta,aTn,bTn,kTn
   real(rk),INTENT(in) :: temp  ! Temperature
!
!LOCALS
   real(rk)  :: fT        !-- Value of the temperature function
   real(rk),PARAMETER  :: tp = 20.0
!
!-------------------------------------------------------------------------------
!BEGIN
   fT = 1.0_rk

   IF ( method /= 1 ) RETURN

   IF (temp > T_max) THEN
       fT = 0.0_rk
   ELSEIF ( temp < T_std ) THEN
       IF (ABS(temp-tp) > 1+MINEXPONENT(temp)/2) THEN
         fT = theta**(temp-tp)
       ENDIF
   ELSE
      IF (ABS(temp-tp) > 1 + MINEXPONENT(temp)/2 .AND. &
          ABS((kTn*(temp-aTn)) + bTn) > 1 + MINEXPONENT(temp)/2) THEN
        fT = theta**(temp-tp) - theta**(kTn*(temp - aTn)) + bTn
       ENDIF
   ENDIF

END FUNCTION fTemp_function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!#                                                                             #
!# A fortran implementation of the quicksort algorithm.                        #
!#                                                                             #
!###############################################################################
RECURSIVE SUBROUTINE qsort(RA,IA,start,end)
!-------------------------------------------------------------------------------
!ARGUMENTS
   real(rk),INTENT(in) :: RA(:)
   INTEGER,INTENT(inout) :: IA(:)
   INTEGER,INTENT(in) :: start,end
!
!LOCALS
  INTEGER :: p, l, r
!
!-------------------------------------------------------------------------------
!BEGIN
  IF ( start .LT. end ) THEN
     l=start+1
     r=end
     p = IA(start);

     DO WHILE(l<r)
        IF (cmp(IA(l), p) .LE. 0) THEN
           l=l+1;
        ELSEIF (cmp(IA(r), p) .GE. 0) THEN
           r=r-1
        ELSE
           CALL swap(IA(l), IA(r))
        ENDIF
     ENDDO
     IF (cmp(IA(l), p) .LT. 0 ) THEN
        CALL swap(IA(l), IA(start))
        l=l-1
     ELSE
        l=l-1
        CALL swap(IA(l), IA(start))
     ENDIF

     CALL qsort(RA,IA,start,l)
     CALL qsort(RA,IA,r,end)
  ENDIF

CONTAINS

   !############################################################################
   SUBROUTINE swap(a, b)
   !----------------------------------------------------------------------------
     INTEGER,intent(inout) :: a, b
     INTEGER t
   !----------------------------------------------------------------------------
   !BEGIN
     t = a
     a = b
     b = t
   END SUBROUTINE swap

   !############################################################################
   INTEGER FUNCTION cmp(l,r)
   !----------------------------------------------------------------------------
      INTEGER,INTENT(in)::l,r
   !----------------------------------------------------------------------------
   !BEGIN
      IF ( RA(l) .LT. RA(r) ) THEN
         cmp = -1
      ELSEIF ( RA(l) .EQ. RA(r) ) THEN
         cmp = 0
      ELSE
         cmp = 1
      ENDIF
   END FUNCTION cmp
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE qsort
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_util
#endif

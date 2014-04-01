!###############################################################################
!#                                                                             #
!# aed_phyto_utils.F90                                                         #
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
!# Created August 2011                                                         #
!#                                                                             #
!###############################################################################

#include "aed.h"

MODULE aed_phyto_utils
!-------------------------------------------------------------------------------
!  aed_phyto_utils --- utility functions for phytoplankton biogeochemical model
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,ONLY : find_free_lun, &
                       exp_integral,  &
                       aed_bio_temp_function, &
                       fTemp_function

   IMPLICIT NONE

!  PRIVATE   ! By default make everything private
!
!  PUBLIC aed_type_phytoplankton
!
   TYPE phyto_data
      ! General Attributes
      CHARACTER(64) :: p_name
      AED_REAL :: p0, Xcc, kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
      ! Growth rate parameters
      INTEGER  :: fT_Method
      AED_REAL :: R_growth, theta_growth, T_std, T_opt, T_max, kTn, aTn, bTn
      ! Light configuration and parameters
      INTEGER  :: lightModel
      AED_REAL :: I_K, I_S, KePHY
      ! Respiration parameters
      AED_REAL :: f_pr, R_resp, k_fdom, k_fres, theta_resp
      ! Salinity parameters
      INTEGER  :: salTol
      AED_REAL :: S_bep, S_maxsp, S_opt
      ! Nitrogen parameters
      INTEGER  :: simDINUptake, simDONUptake, simNFixation, simINDynamics
      AED_REAL :: N_o, K_N, X_nmin, X_nmax, X_ncon, R_nuptake, k_nfix, R_nfix
      ! Phosphorus parameters
      INTEGER  :: simDIPUptake, simIPDynamics
      AED_REAL :: P_0, K_P, X_pmin, X_pmax, X_pcon, R_puptake
      ! Silica parameters
      INTEGER  :: simSiUptake
      AED_REAL :: Si_0, K_Si, X_sicon
      ! Carbon parameters
      INTEGER  :: simCUptake, dic_mode
      ! Sedimentation parameters
      AED_REAL :: w_p
      INTEGER  :: w_model
   END TYPE


   TYPE phyto_nml_data
      CHARACTER(64) :: p_name
      AED_REAL :: p_initial
      AED_REAL :: p0, w_p, Xcc, R_growth !i_min,rmax,alpha,rpn,rpdu,rpdl
      INTEGER  :: fT_Method
      AED_REAL :: theta_growth, T_std, T_opt, T_max
      INTEGER  :: lightModel
      AED_REAL :: I_K, I_S, KePHY
      ! Respiration parameters
      AED_REAL :: f_pr, R_resp, theta_resp, k_fres, k_fdom
      ! Salinity parameters
      INTEGER  :: salTol
      AED_REAL :: S_bep, S_maxsp, S_opt
      ! Nitrogen parameters
      INTEGER  :: simDINUptake, simDONUptake, simNFixation, simINDynamics
      AED_REAL :: N_o, K_N, X_ncon, X_nmin, X_nmax, R_nuptake, k_nfix, R_nfix
      ! Phosphorus parameters
      INTEGER  :: simDIPUptake, simIPDynamics
      AED_REAL :: P_0, K_P, X_pcon, X_pmin, X_pmax, R_puptake
      ! Silica parameters
      INTEGER  :: simSiUptake
      AED_REAL :: Si_0, K_Si, X_sicon
    !  ! Carbon parameters
    !  INTEGER  :: simCUptake, dic_mode
   END TYPE

!Module Locals
   INTEGER  :: ino3, inh4,idon, in2, ifrp, idop

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE phyto_internal_phosphorus(phytos,group,npup,phy,IP,primprod,        &
                                                 fT,pup,respiration,exudation, &
                                                     uptake,excretion,mortality)
!-------------------------------------------------------------------------------
! Calculates the internal phosphorus stores and fluxes
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data),DIMENSION(:),INTENT(in)    :: phytos
   INTEGER,INTENT(in)                          :: group
   INTEGER,INTENT(in)                          :: npup
   AED_REAL,INTENT(in)                         :: phy
   AED_REAL,INTENT(in)                         :: IP
   AED_REAL,INTENT(in)                         :: primprod
   AED_REAL,INTENT(in)                         :: fT,pup,respiration,exudation
   AED_REAL,INTENT(out)                        :: uptake(:),excretion,mortality
!CONSTANTS
   AED_REAL,PARAMETER :: one_e_neg5 = 1e-5
!LOCALS
   AED_REAL :: dumdum1,dumdum2,theX_pcon
   INTEGER  :: c
!
!-------------------------------------------------------------------------------
!BEGIN
   uptake     = zero_
   excretion  = zero_
   mortality  = zero_

   ! Uptake of phosphorus
   IF (phytos(group)%simIPDynamics == 0 .OR. phytos(group)%simIPDynamics == 1) THEN

      ! Static phosphorus uptake function
      ! uptake = X_pcon * mu * phy

      theX_pcon = phytos(group)%X_pcon * phy
      DO c = 1,npup
         ! uptake is spread over relevant sources (assumes evenly)
         uptake(c) = - (theX_pcon/npup) * primprod
      END DO

   ELSEIF (phytos(group)%simIPDynamics == 2) THEN

      ! Dynamic phosphorus uptake function
      ! uptake = R_puptake*fT*phy* (X_pmax-IP/phy)/(X_pmax-X_pmin) * (PO4/(K_P + PO4))

      theX_pcon = IP
      dumdum1  = phytos(group)%R_puptake * fT * phy
      dumdum2  = MAX(one_e_neg5, phytos(group)%X_pmax - (IP / phy))
      dumdum1  = dumdum1 * dumdum2 / (phytos(group)%X_pmax-phytos(group)%X_pmin)
      uptake(1)= -dumdum1 * phyto_fP(phytos,group,frp=pup)
      uptake(2)= 0.0

   ELSE

      ! Unknown phosphorus uptake function
      print *,'STOP: unknown simIPDynamics (',phytos(group)%simIPDynamics,') for: ',phytos(group)%p_name
      STOP

   ENDIF

   ! Release of phosphorus due to excretion from phytoplankton and
   ! contribution of mortality and excretion to OM

   excretion = (respiration*phytos(group)%k_fdom + exudation)*theX_pcon
   mortality = respiration*(1.0-phytos(group)%k_fdom)*theX_pcon


END SUBROUTINE phyto_internal_phosphorus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE phyto_internal_nitrogen(phytos,group,do_N2uptake,phy,IN,primprod,   &
                                   fT,no3up,nh4up,a_nfix,respiration,exudation,&
                                   PNf,uptake,excretion,mortality)
!-------------------------------------------------------------------------------
! Calculates the internal nitrogen stores and fluxes
!-------------------------------------------------------------------------------

!ARGUMENTS
   TYPE(phyto_data),DIMENSION(:),INTENT(in)    :: phytos
   INTEGER,INTENT(in)                          :: group
   LOGICAL,INTENT(in)                          :: do_N2uptake
   AED_REAL,INTENT(in)                         :: phy
   AED_REAL,INTENT(in)                         :: IN
   AED_REAL,INTENT(in)                         :: primprod
   AED_REAL,INTENT(in)                         :: fT,no3up,nh4up
   AED_REAL,INTENT(out)                        :: a_nfix
   AED_REAL,INTENT(in)                         :: respiration,exudation
   AED_REAL,INTENT(out)                        :: PNf
   AED_REAL,INTENT(out)                        :: uptake(:),excretion,mortality
!
!CONSTANTS
   AED_REAL,PARAMETER :: one_e_neg5 = 1e-5
!
!LOCALS
   AED_REAL  :: dumdum1,dumdum2,theX_ncon
!
!-------------------------------------------------------------------------------
!BEGIN
   uptake     = zero_
   excretion  = zero_
   mortality  = zero_


   ! Uptake of nitrogen
   IF (phytos(group)%simINDynamics == 0 .OR. phytos(group)%simINDynamics == 1) THEN

      ! Static nitrogen uptake function (assuming fixed stoichiometry)
      ! uptake = X_ncon * mu * phy

      theX_ncon = phytos(group)%X_ncon * phy
      uptake(1)  = -theX_ncon * primprod

   ELSEIF (phytos(group)%simINDynamics == 2) THEN

      ! Dynamic nitrogen uptake function
      ! uptake = R_nuptake*fT*phy* (X_nmax-IN/phy)/(X_nmax-X_nmin) * (DIN/(K_N + DIN))

      theX_ncon = IN
      dumdum1  = phytos(group)%R_nuptake * fT * phy
      dumdum2  = MAX(phytos(group)%X_nmax - (IN / phy),one_e_neg5)
      dumdum1  = dumdum1 * dumdum2 / (phytos(group)%X_nmax-phytos(group)%X_nmin)
      uptake(1) = dumdum1 * phyto_fN(phytos,group,din=no3up+nh4up)
      uptake(1) = -uptake(1)
   ELSE
      ! Unknown nitrogen uptake function
      print *,'STOP: unknown simINDynamics (',phytos(group)%simINDynamics,') for: ',phytos(group)%p_name
      STOP
   ENDIF

   ! Now find out how much of this was due to N Fixation:
   IF (phytos(group)%simNFixation /= 0) THEN
      a_nfix = phytos(group)%R_nfix * a_nfix * phy
      IF (a_nfix > uptake(1)) THEN
         ! Extreme case:
         a_nfix = uptake(1)
         uptake(1) = zero_
      ELSE
         ! Reduce nuptake by the amount fixed:
         uptake(1) = uptake(1) * (uptake(1)-a_nfix) / uptake(1)
      ENDIF
   ENDIF

   ! Disaggregate N sources to NO3, NH4, DON and N2, based on configuraiton
   PNf = phyto_pN(phytos,group,nh4up,no3up)

   IF (phytos(group)%simDINUptake /= 0) THEN
!     uptake(inh4) = uptake(1) * (1.0-PNf) !inh4 == 2
!     uptake(ino3) = uptake(1) * PNf       !ino3 == 1
      uptake(inh4) = uptake(1) * (PNf)     !inh4 == 2
      uptake(ino3) = uptake(1) * (1.-PNf)  !ino3 == 1
   ENDIF
   IF (phytos(group)%simDONUptake /= 0) THEN
      uptake(idon) = 0.0  !MH to fix  (idon == 3)
   ENDIF
   IF (phytos(group)%simNFixation /= 0 .AND. do_N2uptake) THEN
      uptake(iN2) = a_nfix       ! iN2 == 4
   ENDIF


   ! Release of nitrogen due to excretion from phytoplankton and
   ! contribution of mortality and excretion OM:
   ! (/day +/day)* mg N/ mg C * mgC

   excretion = (respiration*phytos(group)%k_fdom + exudation)*theX_ncon
   mortality = respiration*(1.0-phytos(group)%k_fdom)*theX_ncon

   ! should check here e or m is not exceeding X_nmin

END SUBROUTINE phyto_internal_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_fN(phytos, group, IN, din, don) RESULT(fN)
!-------------------------------------------------------------------------------
! Nitrogen limitation of phytoplankton.
! Michaelis-Menton type formulation or droop model for species with IN
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data),DIMENSION(:),INTENT(in)    :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in),OPTIONAL                :: IN
   AED_REAL,INTENT(in),OPTIONAL                :: din
   AED_REAL,INTENT(in),OPTIONAL                :: don
!
!LOCALS
   AED_REAL :: fN
   AED_REAL :: nup

!-------------------------------------------------------------------------------
!BEGIN
   fN=one_

   IF (PRESENT(din) .OR. PRESENT(don)) THEN
     ! Calculate external nutrient limitation factor
     nup = 0.0
     IF (PRESENT(din) .AND. phytos(group)%simDINUptake == 1) THEN
       nup = nup + din
     ENDIF
     IF (PRESENT(don) .AND. phytos(group)%simDONUptake == 1) THEN
       nup = nup + don
     ENDIF
     fN = (nup-phytos(group)%N_o) / &
           (nup-phytos(group)%N_o+phytos(group)%K_N)
   ELSE
     ! Calculate internal nutrient limitation factor
     fN =   phytos(group)%X_nmax*(1.0-phytos(group)%X_nmin/IN) / &
            (phytos(group)%X_nmax-phytos(group)%X_nmin)
   ENDIF

   IF ( fN < zero_ ) fN = zero_
END FUNCTION phyto_fN
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_fP(phytos, group, IP, frp) RESULT(fP)
!-------------------------------------------------------------------------------
! Phosphorus limitation of phytoplankton
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data),DIMENSION(:),INTENT(in)    :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in), OPTIONAL               :: IP
   AED_REAL,INTENT(in), OPTIONAL               :: frp
!
!LOCALS
   AED_REAL :: fP
!
!-------------------------------------------------------------------------------
!BEGIN
   fP=one_

   IF(PRESENT(frp)) THEN
     fP = (frp-phytos(group)%P_0) / &
             (phytos(group)%K_P + (MAX(zero_, (frp-phytos(group)%P_0))))
   ELSE
     fP = phytos(group)%X_pmax * (1.0 - phytos(group)%X_pmin/IP) / &
                          (phytos(group)%X_pmax-phytos(group)%X_pmin)
   ENDIF

   IF( fP<zero_ ) fP=zero_

END FUNCTION phyto_fP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_fSi(phytos, group, Si) RESULT(fSi)
!-------------------------------------------------------------------------------
! Silica limitation (eg. for diatoms)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data),DIMENSION(:),INTENT(in)    :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in)                         :: Si
!
!LOCALS
   AED_REAL :: fSi
!
!-------------------------------------------------------------------------------
!BEGIN
   fSi = one_

   IF (phytos(group)%simSiUptake == 1) THEN
     fSi = (Si-phytos(group)%Si_0) / &
           (Si-phytos(group)%Si_0+phytos(group)%K_Si)
     IF ( fSi < zero_ ) fSi=zero_
   ENDIF

END FUNCTION phyto_fSi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_pN(phytos,group,NH4,NO3) RESULT(pN)
!-------------------------------------------------------------------------------
! Calculates the relative preference of uptake by phytoplankton of
! ammonia uptake over nitrate.
!-------------------------------------------------------------------------------
   TYPE(phyto_data),DIMENSION(:),INTENT(in)    :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in)                         :: NH4
   AED_REAL,INTENT(in)                         :: NO3
!
!LOCALS
   AED_REAL :: pN
!
!-------------------------------------------------------------------------------
!BEGIN
   pN = zero_

   IF (NH4 > 0.0) THEN
      pN = NH4*NO3 / ((NH4+phytos(group)%K_N)*(NO3+phytos(group)%K_N)) &
         + NH4*phytos(group)%K_N / ((NH4+NO3)*(NO3+phytos(group)%K_N))
   ENDIF
END FUNCTION phyto_pN
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION findMin(a1,a2,a3,a4) RESULT(theMin)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: a1,a2,a3,a4
!LOCALS
   AED_REAL     :: theMin
!
!-------------------------------------------------------------------------------
!BEGIN
   theMin = a1
   IF(a2 < theMin)      theMin = a2
   IF(a3 < theMin)      theMin = a3
   IF(a4 < theMin)      theMin = a4

   IF( theMin<zero_ )  theMin=zero_

END FUNCTION findMin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_respiration(phytos,group,temp) RESULT(respiration)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data),DIMENSION(:),INTENT(in)    :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in)                         :: temp
!
!LOCALS
   AED_REAL :: respiration ! Returns the phytoplankton respiration.
!
!-------------------------------------------------------------------------------
!BEGIN
   respiration = phytos(group)%R_resp * phytos(group)%theta_resp**(temp-20.0)

END FUNCTION phyto_respiration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_salinity(phytos,group,salinity) RESULT(fSal)
!-------------------------------------------------------------------------------
! Salinity tolerance of phytoplankton
! Implmentation based on Griffin et al 2001; Robson and Hamilton, 2004
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data),DIMENSION(:),INTENT(in)    :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in)                         :: salinity
!
!LOCALS
   AED_REAL :: fSal ! Returns the salinity function
   AED_REAL :: tmp1,tmp2,tmp3
   AED_REAL,PARAMETER :: wq_one = 1.0
!
!-------------------------------------------------------------------------------
!BEGIN

   IF (phytos(group)%salTol == 0) THEN
      fSal = 1.0
   ELSEIF (phytos(group)%salTol == 1) THEN
      !# f(S) = 1 at S=S_opt, f(S) = S_bep at S=S_maxsp.
      tmp1 = (phytos(group)%S_bep-1.0) / ((phytos(group)%S_maxsp - phytos(group)%S_opt)**2.0)
      tmp2 = (phytos(group)%S_bep-1.0) * 2.0*phytos(group)%S_opt / &
            ((phytos(group)%S_maxsp-phytos(group)%S_opt)**2.0)
      tmp3 = (phytos(group)%S_bep-1.0) * phytos(group)%S_opt*phytos(group)%S_opt / &
            ((phytos(group)%S_maxsp-phytos(group)%S_opt)**2.0) + 1.0
      IF (salinity>phytos(group)%S_opt) THEN
         fSal = tmp1*(salinity**2.0)-tmp2*salinity+tmp3
      ELSE
         fSal = 1.0
      ENDIF
   ELSEIF (phytos(group)%salTol == 2) THEN
      !# f(S) = 1 at S=S_opt, f(S) = S_bep at S=0.
      IF (salinity<phytos(group)%S_opt) THEN
         fSal = (phytos(group)%S_bep-1.0) * (salinity**2.0)/(phytos(group)%S_opt**2.0) -  &
                      2.0*(phytos(group)%S_bep-1.0)*salinity/phytos(group)%S_opt+phytos(group)%S_bep
      ELSE
        fSal = 1.0
      ENDIF
   ELSEIF (phytos(group)%salTol == 3) THEN
      ! f(S) = 1 at S=S_opt, f(S) = S_bep at S=0 and 2*S_opt.
      IF (salinity < phytos(group)%S_opt) THEN
      fSal = (phytos(group)%S_bep-1.0)*(salinity**2.0)/(phytos(group)%S_opt**2.0)-  &
                      2.0*(phytos(group)%S_bep-1.0)*salinity/phytos(group)%S_opt+phytos(group)%S_bep
      ENDIF
      IF ((salinity>phytos(group)%S_maxsp) .AND. (salinity<(phytos(group)%S_maxsp + phytos(group)%S_opt))) THEN
         fSal = (phytos(group)%S_bep - wq_one)*(phytos(group)%S_maxsp + phytos(group)%S_opt - salinity)**2  &
             / (phytos(group)%S_opt**2) -                                                                             &
             2 * (phytos(group)%S_bep - wq_one) * (phytos(group)%S_maxsp + phytos(group)%S_opt - salinity)  &
             / phytos(group)%S_opt + phytos(group)%S_bep
      ENDIF
      IF ( (salinity >= phytos(group)%S_opt) .AND. (salinity <= phytos(group)%S_maxsp) ) fSal = 1
      IF ( salinity >= (phytos(group)%S_maxsp + phytos(group)%S_opt) ) fSal = phytos(group)%S_bep
   ELSE
      PRINT *,'STOP: Unsupported salTol flag for group: ',group,'=', phytos(group)%salTol
   ENDIF

   IF( fSal < zero_ ) fSal = zero_

END FUNCTION phyto_salinity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_light(phytos, group, par, extc, Io, dz) RESULT(fI)
!-------------------------------------------------------------------------------
! Light limitation of pytoplankton via various model approaches. Refer to
! overview presented in Table 1 of:
!
! Baklouti, M., Diaz, F., Pinazo, C., Faure, V., Quéguiner, B., 2006.
!  Investigation of mechanistic formulations depicting phytoplankton dynamics for
!    models of marine pelagic ecosystems and description of a new model.
!  Progress in Oceanography 71 (1), 1-33.
!
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(phyto_data),DIMENSION(:),INTENT(in)    :: phytos
   INTEGER,INTENT(in)                          :: group
   AED_REAL,INTENT(in)                         :: par
   AED_REAL,INTENT(in)                         :: extc
   AED_REAL,INTENT(in)                         :: Io
   AED_REAL,INTENT(in)                         :: dz
!
!CONSTANTS
   AED_REAL,PARAMETER :: one_e_neg3 = 1e-3
!
!LOCALS
   AED_REAL :: fI !-- Returns the light limitation
   AED_REAL :: par_t,par_b,par_c
   AED_REAL :: z1,z2
   AED_REAL :: x
   AED_REAL, PARAMETER :: A = 5.0, eps = 0.5
!
!-------------------------------------------------------------------------------
!BEGIN
   fI    = 0.0

   ! MH fix this
   par_t = par
   par_b = par_t * EXP( -extc * dz )
   par_c = par_t * EXP( -extc * dz/2. )

   SELECT CASE (phytos(group)%lightModel)
      CASE ( 0 )
         ! Light limitation without photoinhibition.
         ! This is the Webb et al (1974) model solved using the numerical
         ! integration approach as in CAEDYM (Hipsey and Hamilton, 2008)

         IF (Io == zero_) RETURN

         z1 = -par_t / phytos(group)%I_K
         z2 = -par_b / phytos(group)%I_K

         z1 = exp_integral(z1)
         z2 = exp_integral(z2)

         fI = 1.0 + (z2 - z1) / MAX(extc * dz,one_e_neg3)

         ! A simple check
         IF (par_t < 5e-5 .OR. fI < 5e-5) fI = 0.0

      CASE ( 1 )
         ! Light limitation without photoinhibition.
         ! This is the Monod (1950) model.

         x = par_c/phytos(group)%I_K
         fI = x / (one_ + x)

      CASE ( 2 )
         ! Light limitation with photoinhibition.
         ! This is the Steele (1962) model.

         x = par_c/phytos(group)%I_S
         fI = x * EXP(one_ - x)
         IF (par_t < 5e-5 .OR. fI < 5e-5) fI = 0.0

      CASE ( 3 )
         ! Light limitation without photoinhibition.
         ! This is the Webb et al. (1974) model.

         x = par_c/phytos(group)%I_K
         fI = one_ - EXP(-x)

      CASE ( 4 )
         ! Light limitation without photoinhibition.
         ! This is the Jassby and Platt (1976) model.

         x = par_c/phytos(group)%I_K
         fI = TANH(x)

      CASE ( 5 )
         ! Light limitation without photoinhibition.
         ! This is the Chalker (1980) model.

         x = par_c/phytos(group)%I_K
         fI = (EXP(x * (one_ + eps)) - one_) / &
              (EXP(x * (one_ + eps)) + eps)

      CASE ( 6 )
         ! Light limitation with photoinhibition.
         ! This is the Klepper et al. (1988) / Ebenhoh et al. (1997) model.
         x = par_c/phytos(group)%I_S
         fI = ((2.0 + A) * x) / ( one_ + (A * x) + (x * x) )

      CASE ( 7 )
         ! Light limitation with photoinhibition.
         ! This is an integrated form of Steele model.

         fI = ( EXP(1-par_b/phytos(group)%I_S) - &
                EXP(1-par_t/phytos(group)%I_S)   ) / (extc * dz)
   END SELECT

   IF ( fI < zero_ ) fI = zero_
END FUNCTION phyto_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_phyto_utils

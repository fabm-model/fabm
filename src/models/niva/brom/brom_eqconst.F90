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
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev, Jorn Bruggeman
!

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_eqconst
!     Variable identifiers
      type (type_diagnostic_variable_id)   :: id_Kc1,id_Kc2,id_Kw,id_Kb,id_Kp1,id_Kp2,id_Kp3,id_Kc0,id_KSi,id_Knh4,id_Kh2s1,id_Kh2s2
      type (type_dependency_id)            :: id_temp,id_salt
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
   call self%register_diagnostic_variable(self%id_Kh2s1,'Kh2s1','-','[H+][HS-]/[H2S]')
   call self%register_diagnostic_variable(self%id_Kh2s2,'Kh2s2','-','[H+][S2-]/[HS-]')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
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
   real(rk) :: temp,salt
   real(rk) :: Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4,Kh2s1,Kh2s2
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_salt,salt)              ! salinity   

   ! calculate constants needed for the alkalinity components
    call EQCONST(temp,salt,Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4,Kh2s1,Kh2s2)

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
   _SET_DIAGNOSTIC_(self%id_Kh2s1,Kh2s1)
   _SET_DIAGNOSTIC_(self%id_Kh2s2,Kh2s2)
   
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!===========================================================================
!---------------------------------------------------------------------------
!
!               E  Q  C  O  N  S  T       !  for pH-TOTAL !
!
!---------------------------------------------------------------------------
      SUBROUTINE EQCONST(T,S,Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4,Kh2s1,Kh2s2)
 !     call EQCONST(z(k),tem2(k,julianday),sal2(k,julianday),Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4,Kh2s1,Kh2s2)
 
!---------------------------------------------------------------------------
!      default: equilibrium constants after Dickson and Goyet, (1994)
!      "Handbook of methods for the analysis of the various parameters
!      of the carbon dioxide system in sea water", ver. 2, ORNL/CDIAC-74
!      http://cdiac.esd.ornl.gov/oceans/handbook.html.
!       The same values are given in (Grasshoff et al., 1999)
!# if defined co2_const_dickson
!      equilibrium constants Kc1, Kc1 (0>= S <= 40 psu)
!      after Dickson und Millero (1987), Deep Sea Research, Vol. 34
!      corrected constants after Hansson & Mehrbach et. al.
!      check for T=25C, S=35psu: log10(Kc1)=-13.4847, log10(Kc2)=-20.5504,
!# else
!      check for T=25C, S=35psu: ln(Kc1)=-13.4847, ln(Kc2)=-20.5504,
!# endif
!      ln(Kw)=-30.434, ln(Kb)=-19.7964
!# if defined co2_phosphate
!      ln(Kp1)=-3.71, ln(Kp2)=-13.727, ln(Kp2)=-20.24
!# endif
!----------------------------------------------------------------------
      implicit none
      real(rk), intent(in) :: T,S
      real(rk), intent(out) :: Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4,Kh2s1,Kh2s2
      
      real(rk) :: TK,TKR,TLOG,SQ,SLOG,I, Cl
      
! terms used more than once       

      TK   = T + 273.15             ! tk
      TKR  = 1./TK                  ! invtk
      TLOG = LOG(TK)                ! logtk 
      SQ   = SQRT(S)                !sqrts
      SLOG = LOG(1.0e0-0.001005e0*S)
      Cl   = S/1.80655
      I    = 19.924*S/(1000.-1.005*S)  !is

!====================================================================
!    CO2 solubility in water 
!C-------------------------------------------------------------------
!  Kc0 - Henry's constant 
!% Weiss, R. F., Marine Chemistry 2:203-215, 1974.
      Kc0 = EXP(-60.2409+9345.17/TK+23.3585*log((TK)/100.) &
          +S*(0.023517-0.023656*TK/100.+0.0047036*((TK/100.)*(TK/100.))))
!C-------------------------------------------------------------------
     
!====================================================================
     !   CARBONIC ACID
!C-------------------------------------------------------------------
!  Kc1 = [H+][HCO3-]/[H2CO3]
!  (Roy et al., 1993) artificial seawater pH-TOTAL 
      Kc1 = EXP( -2307.1266*TKR + 2.83655 - 1.5529413*TLOG    &
              + (-4.0484*TKR - 0.20760841 - 0.00654208*S)*SQ &
              + 0.08468345*S + SLOG )
!  check for T=25C, S=35psu:  ln(Kc1)= -13.4847   -> CHECKED 091001           
!--------------------------------------------------------------------
!  Kc2 = [H+][CO3--]/[HCO3-]
!  (Roy et al., 1993) artificial seawater pH-TOTAL
      Kc2 = EXP( -3351.6106*TKR - 9.226508 - 0.2005743*TLOG     &
              + (-23.9722*TKR - 0.106901773 - 0.00846934*S)*SQ &
              + 0.1130822*S + SLOG )
!  check for T=25C, S=35psu:  ln(Kc2)= -20.5504   -> CHECKED 091001  
     
! (Dickson)
!      Kc1 = 10**((840.39*TKR - 19.894 + 3.0189*TLOG)*SQ - 0.00668*S &
!               - 6320.81*TKR + 126.3405 - 19.568*TLOG)

!      Kc2 = 10**((690.59*TKR - 17.176 + 2.6719*TLOG)*SQ - 0.0217*S &
!               - 5143.69*TKR + 90.1833 - 14.613*TLOG)

! (Millero p.664 (1995) using Mehrbach et al. data on seawater scale
     !      ak1(i,j,bi,bj)=10.**(-1. _d 0*(3670.7 _d 0*invtk -
     !&          62.008 _d 0 + 9.7944 _d 0*logtk -
     !&          0.0118 _d 0 * s + 0.000116 _d 0*s2))
     !      ak2(i,j,bi,bj)=10.**(-1. _d 0*(1394.7 _d 0*invtk+ 4.777 _d 0-
     !&          0.0184 _d 0*s + 0.000118 _d 0*s2))
!C-------------------------------------------------------------------
     
!====================================================================
!        WATER
!C-------------------------------------------------------------------
!  Kw = [H+][OH-]/H2O
!  DOE(1994) artificial seawater pH-TOTAL
      Kw = EXP( -13847.26*TKR + 148.9652 - 23.6521*TLOG &
              + (118.67*TKR - 5.977 + 1.0495*TLOG)*SQ - 0.01615*S )
!  check for T=25C, S=35psu:  lnKW -30.434       -> CHECKED 091001                           
!C-------------------------------------------------------------------
     
!====================================================================
!        BORIC ACID
!C-------------------------------------------------------------------
!  Kb = [H+][B(OH)4-]/[B(OH)3]      
!  Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
      Kb = EXP( (-8966.90 + (1.728*S - 2890.53)*SQ        &
                          - (77.942 + 0.0996*S)*S)*TKR    &
              + (148.0248 + 137.1942*SQ + 1.62142*S)      &
              - ( 24.4344 +  25.0850*SQ + 0.24740*S)*TLOG &
              + 0.053105*SQ*TK )
!  check for T=25C, S=35psu:  ln(Kb)= -19.7964   -> CHECKED 091001           
!C-------------------------------------------------------------------
     
!====================================================================
!          PHOSPHORIC ACID
!C-------------------------------------------------------------------
!  Kp1 = [H+][H2PO4-]/[H3PO4]
!  DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
      Kp1 = EXP( -4576.752*TKR + 115.525 - 18.453*TLOG &
               + (0.69171 - 106.736*TKR)*SQ - (0.01844 + 0.65643*TKR)*S )
!--------------------------------------------------------------------     
!  Kp2 = [H][HPO4]/[H2PO4]
!  DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
      Kp2 = EXP( -8814.715*TKR + 172.0883 - 27.927*TLOG &
               + (1.35660 - 160.340*TKR)*SQ - (0.05778 - 0.37335*TKR)*S )
!--------------------------------------------------------------------
!  k3p = [H][PO4]/[HPO4]
!  DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)      
      Kp3 = EXP( -3070.75*TKR - 18.141 &
               + (2.81197 + 17.27039*TKR)*SQ - (0.09984 + 44.99486*TKR)*S )
!C-------------------------------------------------------------------
     
!====================================================================
!              SILICIC ACID          
!C-------------------------------------------------------------------
!  KSi = [H+][H3SiO4-]/[Si(OH)4]
!  Millero p.671 (1995) using data from Yao and Millero (1995)
     KSi = EXP( -8904.2/TK+117.385-19.334*TLOG+(-458.79/TK+3.5913)*SQRT(I) &
           + (188.74/TK-1.5998)*I+(-12.1652/TK+0.07871)*I*I + log(1-0.001005*S))
!C-------------------------------------------------------------------

!====================================================================
!             AMMONIA (after Luff et al, 2001)   NH4+ <--Knh4--> NH3 + H+
!C-------------------------------------------------------------------
!C Knh4 = [H+][NH3]/[NH4]  or [NH4][OH-]/[NH3]
!C     note: Millero (1983) gives the DV and DK for pressure correction for the reaction
!C            NH3 + H2O = NH4 + OH
!C     while the value calculated by this program is for the reaction
!C                  NH4 = NH3 + H
! see Luff et al, 2001 how to do it

      Knh4 = EXP(- 6285.33/TK + 0.0001635*TK - 0.25444 &
             + (0.46532 - 123.7184/TK)*SQRT(S) + (-0.01992 + 3.17556/TK)*S)
 !     Knh4=5.2E-6 - proximate value
!C-------------------------------------------------------------------
     
!====================================================================
!C      HYDROGEN SULFIDE (after Luff et al, 2001)  H2S <--Kh2s1--> H+ + HS-
!C                               (Volkov, 1984)    HS- <--Kh2s2--> H+ + S2-
!C-------------------------------------------------------------------
!      LNKS = 225.838 - 13275/TK - 34.6435*LNTK + 0.3449*SQS - 0.0274*S
!      DV = - 11.07 + 0.009*T - 0.942-03*T*T  !<-This is  a corection for pressure! (EY)
!      DK = (-2.89 + 0.054*T)/1000.0
!      LNKS = - DV/(R*TK)*PP + 0.5*DK/(R*TK)*PP*PP + LNKS
!      K(10) = EXP(LNKS)
      Kh2s1= 1000.*EXP(225.838 - 13275./TK - 34.6435*LOG(TK) &
              + 0.3449*SQRT(S) - 0.0274*S)
!     Kh2s1=1.E-7 - proximate value
      Kh2s2=1.E-13 !(Volkov, 1984)
!     Kh2s2=1.E-13 - proximate value
!C-------------------------------------------------------------------
!====================================================================
!              BISULFATE ION          
!C-------------------------------------------------------------------
!C ks = [H+][SO4--]/[HSO4-]
!C dickson (1990, J. chem. Thermodynamics 22, 113)
!           aks(i,j,bi,bj)=exp(-4276.1 _d 0*invtk + 141.328 _d 0 -
!     &          23.093 _d 0*logtk +
!     &   (-13856. _d 0*invtk + 324.57 _d 0 - 47.986 _d 0*logtk)*sqrtis+
!     &   (35474. _d 0*invtk - 771.54 _d 0 + 114.723 _d 0*logtk)*is -
!     &          2698. _d 0*invtk*is**1.5 _d 0 + 1776. _d 0*invtk*is2 +
!     &          log(1.0 _d 0 - 0.001005 _d 0*s))
!C-------------------------------------------------------------------

!====================================================================
!              HYDROGEN FLUORIDE          
!C-------------------------------------------------------------------
!C kf = [H+][F-]/[HF]
!C dickson and Riley (1979) -- change pH scale to total
!           akf(i,j,bi,bj)=exp(1590.2 _d 0*invtk - 12.641 _d 0 +
!     &   1.525 _d 0*sqrtis + log(1.0 _d 0 - 0.001005 _d 0*s) +
!     &   log(1.0 _d 0 + (0.1400 _d 0/96.062 _d 0)*(scl)/aks(i,j,bi,bj)))
     

!OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
! What is needed more...
! Calculate concentrations for borate, sulfate, and fluoride
!C Uppstrom (1974)
!           bt(i,j,bi,bj) = 0.000232 _d 0 * scl/10.811 _d 0
!C Morris & Riley (1966)
!           st(i,j,bi,bj) = 0.14 _d 0 * scl/96.062 _d 0
!C Riley (1965)
!           ft(i,j,bi,bj) = 0.000067 _d 0 * scl/18.9984 _d 0
      return 
      
!-----------------------------------------------------------
      end SUBROUTINE EQCONST
!-----------------------------------------------------------

end module
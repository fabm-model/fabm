#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_carb
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_carb
!     Variable identifiers
      type (type_state_variable_id)        :: id_DIC,id_Alk
      type (type_diagnostic_variable_id)   :: id_pH,id_pCO2,id_Hplus,id_Om_Ca,id_Om_Ar,id_CO3,id_Ca
      type (type_dependency_id)            :: id_temp,id_salt,id_pres
      type (type_dependency_id)            :: id_PO4,id_Si,id_NH4,id_DON,id_H2S,id_Mn3,id_Mn4,id_Fe3
      type (type_dependency_id)            :: id_Kc1,id_Kc2,id_Kw,id_Kb,id_Kp1,id_Kp2,id_Kp3,id_Kc0,id_KSi,id_Knh4,id_Kh2s1,id_Kh2s2
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
   class (type_niva_brom_carb), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%register_state_variable(self%id_DIC, 'DIC', 'mmol/m**3','Total Dissolved Inorganic Carbon', minimum=0.0_rk)
   call self%register_state_variable(self%id_Alk, 'Alk', 'mmol/m**3','Total Alkalinity', minimum=0.0_rk)

   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_pres,standard_variables%pressure)

   call self%register_dependency(self%id_Kc0,'Kc0','-','Henry''s constant')
   call self%register_dependency(self%id_Kc1,'Kc1','-','[H+][HCO3-]/[H2CO3]')
   call self%register_dependency(self%id_Kc2,'Kc2','-','[H+][CO3--]/[HCO3-]')
   call self%register_dependency(self%id_Kw, 'Kw','-','[H+][OH-]/H2O')
   call self%register_dependency(self%id_Kb, 'Kb','-','[H+][B(OH)4-]/[B(OH)3]')
   call self%register_dependency(self%id_Kp1,'Kp1','-','[H+][H2PO4-]/[H3PO4]')
   call self%register_dependency(self%id_Kp2,'Kp2','-','[H][HPO4]/[H2PO4]')
   call self%register_dependency(self%id_Kp3,'Kp3','-','[H][PO4]/[HPO4]')
   call self%register_dependency(self%id_KSi,'KSi','-','[H+][H3SiO4-]/[Si(OH)4]')
   call self%register_dependency(self%id_Knh4,'Knh4','-','[H+][NH3]/[NH4]')
   call self%register_dependency(self%id_Kh2s1,'Kh2s1','-','H2S <--> H+ + HS-')
   call self%register_dependency(self%id_Kh2s2,'Kh2s2','-','HS- <--> H+ + S2-')
   
   call self%register_dependency(self%id_po4,'PO4','mmol/m**3','phosphate')
   call self%register_dependency(self%id_Si, 'Si', 'mmol/m**3','silicate')
   call self%register_dependency(self%id_NH4,'NH4','mmol/m**3','ammonium')
   call self%register_dependency(self%id_DON,'DON','mmol/m**3','dissolved organic nitrogen')
   call self%register_dependency(self%id_H2S,'H2S','mmol/m**3','hydrogen sulfide')
   call self%register_dependency(self%id_Mn3,'Mn3','mmol/m**3','manganese III')
   call self%register_dependency(self%id_Mn4,'Mn4','mmol/m**3','manganese IV')
   call self%register_dependency(self%id_Fe3,'Fe3','mmol/m**3','iron III')
   
   call self%register_diagnostic_variable(self%id_pH,'pH','-','pH')
   call self%register_diagnostic_variable(self%id_pCO2,'pCO2','ppm','partial pressure of CO2')
   call self%register_diagnostic_variable(self%id_Hplus, 'Hplus', 'mmol/m**3','H+ Hydrogen')
   call self%register_diagnostic_variable(self%id_Om_Ca,'Om_Ca','-','CaCO3-Calcite saturation')
   call self%register_diagnostic_variable(self%id_Om_Ar,'Om_Ar','-','CaCO3-Aragonite saturation')
   call self%register_diagnostic_variable(self%id_CO3,'CO3','mmol/m**3','CO3--')
   call self%register_diagnostic_variable(self%id_Ca,'Ca','mmol/m**3','Ca++')

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
   class (type_niva_brom_carb),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:
   real(rk) :: temp,salt,pres
   real(rk) :: DIC,Alk,PO4,Si,NH4,DON,H2S,Mn3,Mn4,Fe3
   real(rk) :: Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4,Kh2s1,Kh2s2
   
   real(rk) :: H_,pH,Om_Ca,Om_Ar
   real(rk) :: co2,pCO2,hco3,co3,Ca
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_salt,salt)              ! salinity   
   _GET_(self%id_pres,pres)              ! pressure in dbar

   ! Our own state variables
   _GET_(self%id_Alk,Alk)
   _GET_(self%id_DIC,DIC)

   ! External (state) variables.
   _GET_(self%id_PO4,PO4)
   _GET_(self%id_Si,Si)
   _GET_(self%id_NH4,NH4)
   _GET_(self%id_DON,DON)
   _GET_(self%id_H2S,H2S)
   _GET_(self%id_Mn3,Mn3)
   _GET_(self%id_Mn4,Mn4)
   _GET_(self%id_Fe3,Fe3)

   ! Equilibrium constants
   _GET_(self%id_Kc1,  Kc1)
   _GET_(self%id_Kc2,  Kc2)
   _GET_(self%id_Kw,   Kw)
   _GET_(self%id_Kb,   Kb)
   _GET_(self%id_Kp1,  Kp1)
   _GET_(self%id_Kp2,  Kp2)
   _GET_(self%id_Kp3,  Kp3)
   _GET_(self%id_Kc0,  Kc0)
   _GET_(self%id_KSi,  KSi)
   _GET_(self%id_Knh4, Knh4)
   _GET_(self%id_Kh2s1,Kh2s1)
   _GET_(self%id_Kh2s2,Kh2s2)

! calculate pH(tot) as a function of TIC and total ALK
      call pHiter (salt,Alk,DIC, &
             PO4,Si,NH4,DON, &
             H2S,Mn3,Mn4,Fe3, &
             Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4,Kh2s1,Kh2s2, &
             H_,pH)             
      
! calculate all the others as a function of pH1(k), alk1(k), tic1(k) and constants
      call CARFIN(temp,salt,0.1_rk*pres,Kc0,Kc1,Kc2, &
             DIC,H_, &
             co2,pCO2,hco3,co3, &
             Ca,Om_Ca,Om_Ar) 
   
   _SET_DIAGNOSTIC_(self%id_pH,pH)
   _SET_DIAGNOSTIC_(self%id_pCO2,pCO2)
   _SET_DIAGNOSTIC_(self%id_Hplus,H_)
   _SET_DIAGNOSTIC_(self%id_Om_Ca,Om_Ca)
   _SET_DIAGNOSTIC_(self%id_Om_Ar,Om_Ar)
   _SET_DIAGNOSTIC_(self%id_CO3,co3)   
   _SET_DIAGNOSTIC_(self%id_Ca,Ca)   
   
! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!============================================================ 
!-----------------------------------------------------------
!
!      p H i t e r     !  pH-TOTAL =f(At,Ct,S(Bt),Pt,Sit,NHt,H2St,Mn3,Mn4,Fe3,DOC)!  EYA 18.02.2014
!
!-----------------------------------------------------------

  subroutine pHiter (S,Alk,Ct,Pt,Sit_,NHt,DOC,H2St,Mn3,Mn4,Fe3,  &
             Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4,Kh2s1,Kh2s2,H_,ph )
!-----------------------------------------------------------
 implicit none
!
 real(rk),intent(in) :: S,Alk,Ct,Pt,Sit_,NHt,DOC,H2St,Mn3,Mn4,Fe3
 real(rk),intent(in) :: Kc1,Kc2,Kw,Kb,Kp1,Kp2,Kp3,Kc0,KSi,Knh4,Kh2s1,Kh2s2
 real(rk), intent(out) :: H_,ph

 real(rk) :: Bt,T1,T2,T12,K12p,K123p,HKR123P,AH,ADH,DH,TK,T,Sit
 integer  iter, oldcolor
!=======================================================================
! TOTAL BORON after Uppstrom, 1974, Deep sea Res., 21, p. 161ff.
! see Dickson, Goyet, Handbook, SOP 3
        Bt  = S * 1.188e-5 ! 1.188e-5*35 = 4.16e-4 
        Sit = 0. ! Sit_
 		!alk1(k)=	Cc(k,i_Alk) &                    ! Total alkalinity     
  !                  - Cc(k,i_DON)*0.10 &  ! OM alkalinity 
  !                  - Cc(k,i_NH4)*Knh4/(Knh4+H_) &   ! ammonia alkalinity 
  !  !                - Cc(k,i_NO2) &
  !                  - Cc(k,i_PO4)*((K12p-H_*H_)*H_+2e0*K123p)*HKR123p & ! Phosphate alkalinity    
  !                  - Cc(k,i_Si)*KSi/(KSi+H_) & ! Silicic alkalinity ! ksilocal / (ksilocal + hguess)
  !                  + Cc(k,i_Mn3) &  ! Mn(III) alkalinity    
  !                  + 2.*Cc(k,i_Mn4) &  ! Mn(IV) alkalinity                      
  !                  + Cc(k,i_Fe3) &  ! Fe(III) alkalinity  
  !                  - Cc(k,i_H2S)*Kh2s1/(Kh2s1+H_) &  ! HS alkalinity 
  !                  - 2.*(Cc(k,i_H2S)*Kh2s1/(Kh2s1+H_))*Kh2s2/(Kh2s2+H_) !S2 alkalinity        
  
!-----------------------------------------------------------  
!-----NEWTON-RAPHSON METHOD
!-----------------------------------------------------------
!     iteration for the solution of pH for given 
!     total alkalinity and total carbon.
!     for example:  x2 = x1 - (f(x1)/f'(x1)) 
! 
!     Ct = [H2CO3] + [HCO3-] + [CO3--] 
!     [HCO3-] = [CO3--]*[H]/Kc2 
!     [H2CO3] = [HCO3-]*[H]/Kc1 = [CO3--]*[H]*[H]/Kc2*Kc1 
!
!     that is 
!     Ct = [CO3--*(1 + [H]/Kc2 +[H]*[H]/Kc1*Kc2)] 
! = > [CO3--] = Ct/(1 + [H]/Kc2 +[H]*[H]/Kc1*Kc2)
! = > [HCO3--] = ([H]/Kc2)*(Ct/(1 + [H]/Kc2 +[H]*[H]/Kc1*Kc2))
! = >   Alk_c = [HCO3--] + 2*[CO3--] = 
! = >   =([H]/Kc2)*(Ct/(1 + [H]/Kc2 +[H]*[H]/Kc1*Kc2)) +
! = >     + 2.*(Ct/(1 + [H]/Kc2 +[H]*[H]/Kc1*Kc2)) =
! = >   = ([H]/Kc2+2.)*Ct/(1 + [H]/Kc2 +[H]*[H]/Kc1*Kc2) 


!     with the estimate start value for [H] get an estimated value for
!     alkalinity A(H).  With the deviation of Alk ( A(H) - Alk ) get a 
!     correction for [H] and than go on with iteration until A(H)=Alk.   
! 
!     Alk = [HCO3-] + 2[CO3--] + [B(OH)4-] + [OH-] - [H+] +- ...
! 
!     KW = [H][OH-]
!     KB = [H][B(OH)4-]/[B(OH)3], Bt = [B(OH)3] + [B(OH)4-]
! 
!     Build the difference between estimated total alkalinity and the 
!     prescribed value for total alkanlinity
!
! run [H+] iteration as outer loop to allow vectorisation of inner loop
! 
!   H_ = 1.e-8      ! first guess
   H_ = 0.5e-7     ! first guess
!C-------------------------------------------------------------------

!C-------------------------------------------------------------------
   do iter=1,100
    
! write (*,*) iter,H_,Ct,Alk,ph
!c set some common combinations of parameters used in
!c the iterative [H+] solvers

      T1  = H_/Kc1
      T2  = H_/Kc2
      T12 = (1e0+T2+T1*T2)
!# if defined co2_phosphate
      K12p  = Kp1*Kp2
      K123p = Kp1*Kp2*Kp3
      HKR123p = 1e0/(((H_+Kp1)*H_+K12p)*H_+K123p)
      
      !AH = Ct*(2.+T2)/T12     &
      !     + Bt*Kb/(Kb+H_)    &       
      !     + Kw/H_ - H_ - Alk &
      !     + Pt*((K12p-H_*H_)*H_+2e0*K123p)*HKR123p     
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW in normal notation:
!   "Alk" =  [HCO3-] + 2[CO3--]  i.e.= ([H]/Kc2+2.)*Ct/(1 + [H]/Kc2 +[H]*[H]/(Kc1*Kc2)) !Carbonate alkalinity
     AH = Ct*(2.+H_/Kc2)/(1.+H_/Kc2+H_/Kc1*H_/Kc2)     &
!           [B(OH)4-]       i.e.= Btot*Kb/(Kb+[H+])       ! boric alkalinity     
           + Bt*Kb/(Kb+H_)    &
!           [OH-]   i.e.= Kw/[H+]
           + Kw/H_ &
!           [H+]
           - H_ &
!           Alk_tot 
           - Alk &
!            [H2PO4-] +  2.*[HPO4--] + 3.*[PO4---]  i.e.  ! phosphoric alkalinity          
           + Pt*((Kp1*Kp2-H_*H_)*H_+2e0*Kp1*Kp2*Kp3) &
                 /(((H_+Kp1)*H_+Kp1*Kp2)*H_+Kp1*Kp2*Kp3)  &
!            [H3SiO4-] i.e.=Sit*KSi/(KSi+[H+])! silicate alkalinity
           + Sit*KSi/(KSi+H_) &      
!            [HS-] i.e.=[H2St] *Kh2s1/(Kh2s1+[H+])    ! hydrogen sulfidic alkalinity          
           + H2St*Kh2s1/(Kh2s1+H_) &        
!            [NH3] i.e.  = NHt*Knh4/(Knh4+[H+])
           + NHt*Knh4/(Knh4+H_)  !&     !ammonia alkalinity
    
!     ADH is the derivative of AH i.e. d(AH)/d[H+]
!---------------------------------------------
      ADH = Ct*(T2-4e0*T1-T12)/(Kc2*T12*T12)    &
            - Bt*Kb/(Kb+H_)**2.          &
            - Kw/H_**2. - 1.                      &
            - Pt*((((Kp1*H_+4e0*K12p)*H_ &
            - Sit*KSi/(KSi+H_)**2. &
            + Kp1*K12p+9e0*K123p)*H_ &
             + 4e0*Kp1*K123p)*H_+K12p*K123p)*HKR123p*HKR123p &    
            - H2St*Kh2s1/(Kh2s1+H_)**2.   &
            - NHt*Knh4/(Knh4+H_)**2.   
!# else
!      ADH = Ct*(T2-4e0*T1-T12)/(Kc2*T12*T12)    &
!            - Bt*Kb/(Kb+H_)**2          &
!            - Kw/H_**2 - 1e0
!# endif
      dh = 0.1*AH/ADH
      H_ = H_ - dh   
      !if (dh<0.000000000001)  then
      !    goto 1 
      !endif 
      !
       !oldcolor = SETTEXTCOLORRGB(Z'000000FF')
       ! ph = -LOG10(H_)
       ! write(*,*) iter,pH
       ! 
   enddo !iter

1      ph = -LOG10(H_)
!pause 123
        return 
!-----------------------------------------------------------
             END SUBROUTINE pHiter
!-----------------------------------------------------------

!===========================================================================
!---------------------------------------------------------------------------
!
!           C  a  C  O  3  s  o  l  u  b
!
!---------------------------------------------------------------------------
 subroutine CaCO3solub(temp, salt, Pbar, Ca, KCal, KAra)
!---------------------------------------------------------------------------
!  input: temp, salt, Pbar
! output: Ca, KCal, KAra
    implicit none
    real(rk),intent(in)  :: temp, salt, Pbar
    real(rk),intent(out) :: Ca, KCal, KAra
    real(rk) :: tempK, logKCal, logKAra, RT, &
        deltaVKCal,KappaKCal,lnKCalfac,deltaVKAra,  &
        KappaKAra, lnKArafac

!% '***********************************************************************
!% ' SUB CaSolubility, version 01.05, 05-23-97, written by Ernie Lewis.
!% ' Inputs: WhichKs%, Sal, TempCi, Pdbari, TCi, pHi, Kc1, Kc2
!% ' Outputs: OmegaCa, OmegaAr
!% ' This calculates omega, the solubility ratio, for calcite and aragonite.
!% ' This is defined by: Omega = [CO3--]*[Ca++]./Ksp,
!% '       where Ksp is the solubility product (either KCa or KAr).
!% '***********************************************************************
!% ' These are from:
!% ' Mucci, Alphonso, The solubility of calcite and aragonite in seawater
!% '       at various salinities, temperatures, and one atmosphere total
!% '       pressure, American Journal of Science 283:781-799, 1983.
!% ' Ingle, S. E., Solubility of calcite in the ocean,
!% '       Marine Chemistry 3:301-319, 1975,
!% ' Millero, Frank, The thermodynamics of the carbonate system in seawater,
!% '       Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
!% ' Ingle et al, The solubility of calcite in seawater at atmospheric pressure
!% '       and 35%o salinity, Marine Chemistry 1:295-307, 1973.
!% ' Berner, R. A., The solubility of calcite and aragonite in seawater in
!% '       atmospheric pressure and 34.5%o salinity, American Journal of
!% '       Science 276:713-730, 1976.
!% ' Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
!% ' Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
!% '       boric acid, and the pHi of seawater, Limnology and Oceanography
!% '       13:403-417, 1968.
!% '***********************************************************************

!F=(WhichKs~=6 & WhichKs~=7);
!if any(F)
!    % CalculateCa:
!    % '       Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
!    % '       this is .010285.*Sali./35
    Ca      = 0.02128/40.087*(salt/1.80655)    ! in mol/kg-SW
    tempK   = temp + 273.15
      
 
!     CalciteSolubility:
!           Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    logKCal =     -171.9065 - 0.077993*tempK + 2839.319/tempK &
                 + 71.595* LOG10(tempK) &
                 + (-0.77712 + 0.0028426*TempK + 178.34/tempK)*salt**(0.5) &
                 - 0.07711*salt + 0.0041249*salt**(1.5)
!  check for T=25C, S=35psu: logKcal = 6.3693
    KCal = 10.**(logKCal)    ! in (mol/kg-SW)^2
! AragoniteSolubility:
!         Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    logKAra =    -171.945 - 0.077993*tempK + 2903.293/tempK &
                 + 71.595* LOG10(tempK) &
                 + (-0.068393 + 0.0017276*tempK + 88.135/tempK)*salt**(0.5) &
                 - 0.10018*salt + 0.0059415*salt**(1.5)
!  check for T=25C, S=35psu: logKcal = 6.1883
    KAra = 10.**(logKAra)    ! in (mol/kg-SW)^2
     !write (*,*) Ca, Tempk,logKCal,KCal,logKAra,KAra
     ! pause 201
!    % PressureCorrectionForCalcite:
!    % '       Ingle, Marine Chemistry 3:301-319, 1975
!    % '       same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
!    % '       has typos (-.5304, -.3692, and 10^3 for Kappa factor)
!    deltaVKCal = -48.76 + 0.5304*temp;
!    KappaKCal  = (-11.76 + 0.3692*temp)/1000;
!    lnKCalfac  = (-deltaVKCal + 0.5*KappaKCal*Pbar)*Pbar/RT;
!    KCal       = KCal*exp(lnKCalfac);

!    % PressureCorrectionForAragonite:
!    % '       Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
!    % '       same as Millero, GCA 1995 except for typos (-.5304, -.3692,
!    % '       and 10^3 for Kappa factor)
!    deltaVKAra = deltaVKCal + 2.8;
!    KappaKAra  = KappaKCal;
!    lnKArafac  = (-deltaVKAra + 0.5*KappaKAra*Pbar)*Pbar/RT;
!    KAra       = KAr*exp(lnKArafac);

        return 
!-----------------------------------------------------------
        END SUBROUTINE CaCO3solub
!-----------------------------------------------------------
             
!=========================================================================== 
!---------------------------------------------------------------------------
!
!               C A R F I N       !  for pH-TOTAL !  EYA 2010-08-17
!
!---------------------------------------------------------------------------
 subroutine CARFIN (temp, salt, Pbar,Kc0,Kc1,Kc2,tic,H_,co2,pco2,hco3,co3,Ca,KsatCal,KsatAra)
!---------------------------------------------------------------------------
! input: temp, salt, Pbar,Kc0,Kc1,Kc2,tic,H_,
! output: co2,pco2,hco3,co3,Ca,KsatCal,KsatAra
 
implicit none
        real(rk),intent(in)  :: temp,salt,Pbar,Kc0,Kc1,Kc2,tic,H_
        real(rk),intent(out) :: co2,pco2,hco3,co3
        real(rk),intent(out) :: Ca,KsatCal,KsatAra
        real(rk) K_Cal,K_Ara
     
           hco3 = tic/(1.+H_/Kc1+Kc2/H_)       !these are in [uM]
           co3  = tic/(1.+H_/Kc2+H_*H_/Kc1/Kc2)
           co2  = tic/(1.+Kc1/H_+Kc1*Kc2/H_/H_)  
           pco2 = co2/Kc0               !  [uatm]  
           call CaCO3solub(temp, salt, Pbar, Ca, K_Cal, K_Ara)
           KsatCal=(co3/1000000.)*Ca/K_Cal  !000.  !Saturation (Omega) for calcite 
           KsatAra=(co3/1000000.)*Ca/K_Ara  !000.  !Saturation (Omega) for aragonite 
       return 
!-----------------------------------------------------------
 END SUBROUTINE CARFIN

end module
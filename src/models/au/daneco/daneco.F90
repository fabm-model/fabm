#ifdef _FABM_F2003_
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_au_daneco Danish Ecologocal model
! Aarhus university, department of BioScience
!
! !INTERFACE:
   module fabm_au_daneco
!
! !DESCRIPTION:
!
!
! !USES:
   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:
!  Original author(s): Janus Larsen
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model), public :: type_au_daneco
!     Variable identifiers
      type (type_state_variable_id)        :: id_B, id_N, id_P !Micro plankton carbon, nitrogen and phosphorous
      type (type_state_variable_id)        :: id_C, id_M, id_W !Detritus carbon, nitrogen and phosphorous
      type (type_state_variable_id)        :: id_NO3, id_NH4, id_PO4 !Nutrients
      type (type_state_variable_id)        :: id_O2 !Oxygen

      type (type_dependency_id)            :: id_par, id_temp, id_salt
      type (type_horizontal_dependency_id) :: id_I_0

      type (type_diagnostic_variable_id)   :: id_dPAR, id_dX

      type (type_conserved_quantity_id)    :: id_totN, id_totP !total nitrogen and phosphorous

!     Model parameters
      real(rk) :: B_initial,N_initial,P_initial,C_initial,M_initial,W_initial,NO3_initial,NH4_initial,PO4_initial,O2_initial
      real(rk) :: G,gam,e,qOB,qONO,qONH,qOC,qh,eta,KNH4,KNO3,KPO4,maxQNB,maxQPB,NOupmax,NHupmax,POupmax
      real(rk) :: qT,Tref,mymax,minQNB,minQPB,qXN,alfa,b0,r0,Crmax,Mrmax,Wrmax,NHrmax,KNHO,KCO,minQMC,minQWC
      real(rk) :: w_d,kw,eS,eC,eX

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
      procedure :: get_light_extinction
      procedure :: get_conserved_quantities
   end type type_au_daneco
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the daneco namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_au_daneco), intent(inout), target :: self
   integer,                intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Karen Timmermann, Janus Larsen
!
! !LOCAL VARIABLES:
   real(rk) :: B_initial
   real(rk) :: N_initial
   real(rk) :: P_initial
   real(rk) :: C_initial
   real(rk) :: M_initial
   real(rk) :: W_initial
   real(rk) :: NO3_initial
   real(rk) :: NH4_initial
   real(rk) :: PO4_initial
   real(rk) :: O2_initial

   real(rk) :: G
   real(rk) :: gam
   real(rk) :: e
   real(rk) :: qOB
   real(rk) :: qONO
   real(rk) :: qONH
   real(rk) :: qOC
   real(rk) :: qh
   real(rk) :: eta
   real(rk) :: KNH4
   real(rk) :: KNO3
   real(rk) :: KPO4
   real(rk) :: maxQNB
   real(rk) :: maxQPB
   real(rk) :: NOupmax
   real(rk) :: NHupmax
   real(rk) :: POupmax
   real(rk) :: qT
   real(rk) :: Tref
   real(rk) :: mymax
   real(rk) :: minQNB
   real(rk) :: minQPB
   real(rk) :: qXN
   real(rk) :: alfa
   real(rk) :: b0
   real(rk) :: r0
   real(rk) :: Crmax
   real(rk) :: Mrmax
   real(rk) :: Wrmax
   real(rk) :: NHrmax
   real(rk) :: KNHO
   real(rk) :: KCO
   real(rk) :: minQMC
   real(rk) :: minQWC

   real(rk) :: w_d !detritus settling rate (m/s)

   real(rk) :: kw  !2.06 !light extinction coef
   real(rk) :: eS  !light extinction coef (salt)
   real(rk) :: eC  !light extinction coef (carbon)
   real(rk) :: eX  !light extinction coef (chlorophyll)

   namelist /au_daneco/ B_initial,N_initial,P_initial,C_initial,M_initial,W_initial,NO3_initial,NH4_initial,PO4_initial,O2_initial, &
                        G,gam,e,qOB,qONO,qONH,qOC,qh,eta,KNH4,KNO3,KPO4,maxQNB,maxQPB,NOupmax,NHupmax,POupmax, &
                        qT,Tref,mymax,minQNB,minQPB,qXN,alfa,b0,r0,Crmax,Mrmax,Wrmax,NHrmax,KNHO,KCO,minQMC,minQWC, &
                        w_d,kw,eS,eC,eX
!EOP
!-----------------------------------------------------------------------
!BOC
   B_initial=1_rk
   N_initial=1_rk
   P_initial=1_rk
   C_initial=1_rk
   M_initial=1_rk
   W_initial=1_rk
   NO3_initial=1_rk
   NH4_initial=1_rk
   PO4_initial=1_rk
   O2_initial=1_rk

   G=0.0_rk !5.79e-7_rk
   gam=0.8_rk
   e=0.5_rk
   qOB=1_rk
   qONO=2_rk
   qONH=2_rk
   qOC=1_rk
   qh=0.18_rk
   eta=0.3_rk
   KNH4=0.24_rk
   KNO3=0.32_rk
   KPO4=0.02_rk
   maxQNB=0.2_rk
   maxQPB=0.013_rk
   NOupmax=4.05e-6_rk
   NHupmax=1.22e-5_rk
   POupmax=2.34e-5_rk
   qT=0.07_rk
   Tref=20_rk
   mymax=2.70e-5_rk
   minQNB=0.02_rk
   minQPB=0.001_rk
   qXN=0.7_rk
   alfa=0.2295_rk
   b0=0.725_rk
   r0=5.1e-7_rk
   Crmax=6.94e-7_rk
   Mrmax=9.26e-7_rk
   Wrmax=9.26e-7_rk
   NHrmax=1.16e-6_rk
   KNHO=30_rk
   KCO=10_rk
   minQMC=0.09_rk
   minQWC=0.005_rk

   w_d=0.0001_rk !detritus settling rate (m/s)

   kw=0.2_rk     !2.06 !light extinction coef
   eS=0.06_rk    !light extinction coef (salt)
   eC=0.002_rk   !light extinction coef (carbon)
   eX=0.03_rk    !light extinction coef (chlorophyll)

   ! Read the namelist
   read(configunit,nml=au_daneco,err=99,end=100)

!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   self%B_initial=B_initial
   self%N_initial=N_initial
   self%P_initial=P_initial
   self%C_initial=C_initial
   self%M_initial=M_initial
   self%W_initial=W_initial
   self%NO3_initial=NO3_initial
   self%NH4_initial=NH4_initial
   self%PO4_initial=PO4_initial
   self%O2_initial=O2_initial
   self%G=G
   self%gam=gam
   self%e=e
   self%qOB=qOB
   self%qONO=qONO
   self%qONH=qONH
   self%qOC=qOC
   self%qh=qh
   self%eta=eta
   self%KNH4=KNH4
   self%KNO3=KNO3
   self%KPO4=KPO4
   self%maxQNB=maxQNB
   self%maxQPB=maxQPB
   self%NOupmax=NOupmax
   self%NHupmax=NHupmax
   self%POupmax=POupmax
   self%qT=qT
   self%Tref=Tref
   self%mymax=mymax
   self%minQNB=minQNB
   self%minQPB=minQPB
   self%qXN=qXN
   self%alfa=alfa
   self%b0=b0
   self%r0=r0
   self%Crmax=Crmax
   self%Mrmax=Mrmax
   self%Wrmax=Wrmax
   self%NHrmax=NHrmax
   self%KNHO=KNHO
   self%KCO=KCO
   self%minQMC=minQMC
   self%minQWC=minQWC

   self%w_d=w_d

   self%kw=kw
   self%eS=eS
   self%eC=eC
   self%eX=eX

!  Register state variables
   call self%register_state_variable(self%id_B,'B','mmol/m**3','Plankton Carbon',B_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_N,'N','mmol/m**3','Plankton Nitrogen',N_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_P,'P','mmol/m**3','Plankton Phophorous',P_initial,minimum=_ZERO_,no_river_dilution=.true.)

   call self%register_state_variable(self%id_C,'C','mmol/m**3','Detritus Carbon',C_initial,minimum=_ZERO_,vertical_movement=w_d)
   call self%register_state_variable(self%id_M,'M','mmol/m**3','Detritus Nitrogen',M_initial,minimum=_ZERO_,vertical_movement=w_d)
   call self%register_state_variable(self%id_W,'W','mmol/m**3','Detritus Phophorous',W_initial,minimum=_ZERO_,vertical_movement=w_d)

   call self%register_state_variable(self%id_NO3,'NO3','mmol/m**3','Nitrate',NO3_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_NH4,'NH4','mmol/m**3','Ammonium',NH4_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_PO4,'PO4','mmol/m**3','Phosphate',PO4_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_O2,'O2','mmol/m**3','Oxygen',O2_initial,minimum=_ZERO_,no_river_dilution=.true.)

!  Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2', 'photosynthetically active radiation',time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_dX,'X','mg/m**3', 'chl.a concentration',time_treatment=time_treatment_averaged)

!  Register conserved quantities
   call self%register_conserved_quantity(self%id_totN,'totN','mmol/m**3','totN')
   call self%register_conserved_quantity(self%id_totP,'totP','mmol/m**3','totP')

!  Register environmental dependencies
   call self%register_dependency(self%id_par,varname_par)
   call self%register_dependency(self%id_I_0,varname_par_sf)
   call self%register_dependency(self%id_temp,varname_temp)
   call self%register_dependency(self%id_salt,varname_salt)

   return

99 call fatal_error('au_daneco_create','Error reading namelist au_daneco.')

100 call fatal_error('au_daneco_create','Namelist au_daneco was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   subroutine do(self,_FABM_ARGS_DO_RHS_)
!
! !INPUT PARAMETERS:
   class (type_au_daneco),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !LOCAL VARIABLES:
   real(rk) :: B,N,P,C,M,W,NO3,NH4,PO4,O2 !local state variabels
   real(rk) :: NOup !Nitrate uptake
   real(rk) :: NHup !Ammonium uptake
   real(rk) :: POup !Phosphorous uptake
   real(rk) :: my   !Growth
   real(rk) :: Cr   !Carbon remineralisation
   real(rk) :: Mr   !Nitrogen remineralisation
   real(rk) :: Wr   !Phosphorous remineralisation
   real(rk) :: NHr  !Ammonium remineralisation

   real(rk) :: par,I_0,temp !local denpendency variables
   real(rk) :: my1  !Growth light
   real(rk) :: my2  !Growth nitrogen
   real(rk) :: my3  !Growth phosphorous
   real(rk) :: QNB  !Plankton nitrogen to carbon ratio
   real(rk) :: QPB  !Plankton phosphorous to carbon ratio
   real(rk) :: fT   !Temperature function
   real(rk) :: X    !Chlorophyll
   real(rk) :: QMC  !Detritus nitrogen to carbon ratio
   real(rk) :: QWC  !Detritus phosphorous to carbon ratio
!#define _USE_DT_
#if _USE_DT_
   real(rk) :: dt
#endif
   real(rk) :: dO2 !change in O2
!EOP
!-----------------------------------------------------------------------
!BOC
!  Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

!  Retrieve current (local) state variable values.
   _GET_(self%id_B,B)
   _GET_(self%id_N,N)
   _GET_(self%id_P,P)
   _GET_(self%id_C,C)
   _GET_(self%id_M,M)
   _GET_(self%id_W,W)
   _GET_(self%id_NO3,NO3)
   _GET_(self%id_NH4,NH4)
   _GET_(self%id_PO4,PO4)
   _GET_(self%id_O2,O2)

!  Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)           ! local water temperature
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation

#if _USE_DT_
   dt=80._rk ! valid for Skivefiord setup only
#endif

   QWC=W/C
   QMC=M/C
   fT=exp(self%qT*(temp-self%Tref))
   if (QMC .ge. self%minQMC) then
      Cr=fT*(O2/(self%KCO+O2))*self%Crmax*(1-self%minQMC/QMC)**2
   else
      Cr=0.0_rk
   end if
   if (QMC .ge. self%minQMC) then
      Mr=fT*self%Mrmax*(1-self%minQMC/QMC)**2
   else
      Mr=0.0_rk
   end if
   if (QWC .ge. self%minQWC) then
      Wr=fT*self%Wrmax*(1-self%minQWC/QWC)**2
   else
      Wr=0.0_rk
   end if
!  To avoid neg. NH4 conc. Limitation on nitrification
   NHr=fT*self%NHrmax*(O2/(self%KNHO+O2))

#if _USE_DT_
   if (NH4.lt.NHr*NH4*dt) then
      NHr = NH4/(NH4*dt)
   end if
#endif

   if (B .gt. 0._rk) then
      QNB=max(N/B,self%minQNB)
   else
      QNB=self%minQNB
   end if

   if (B .gt. 0._rk) then
      QPB=max(P/B,self%minQPB)
   else
      QPB=self%minQPB
   end if

   if (QNB .ge. self%minQNB) then
      my3=self%mymax*fT*(1._rk-self%minQNB/QNB)
   else
      my3=0.0_rk
   end if
   if (QPB .ge. self%minQPB) then
      my2=self%mymax*fT*(1._rk-self%minQPB/QPB)
   else
      my2=0.0_rk
   end if
   X=self%qXN*N
   my1=(self%alfa*par*X-self%r0)/(1._rk+self%b0)
   my=min(my1,min(my2,my3))
   if (my .lt. 0._rk) then
      my1=(self%alfa*par*X-self%r0)
      my=min(my1,min(my2,my3))
   end if

   if (QNB .le. self%maxQNB) then
      NOup=self%NOupmax*fT*NO3/(self%KNO3+NO3)*NH4/(self%KNH4+NH4)*(1._rk-(QNB-self%qh*self%eta)/(self%maxQNB-self%qh*self%eta))
   else
      NOup=0.0_rk
   end if

#if _USE_DT_
   !NOup=if(NO3.ge.NOup*B*dt,NOup,NO3/(B*dt))
   if (NO3.lt.NOup*B*dt) then
      NOup = NO3/(B*dt)
   end if
#endif

   if (QNB .le. self%maxQNB) then
      NHup=self%NHupmax*fT*NH4/(self%KNH4+NH4)*(1._rk-(QNB-self%qh*self%eta)/(self%maxQNB-self%qh*self%eta))
   else
      NHup=0.0_rk
   end if

#if _USE_DT_
   !NHup=if(NH4.ge.NHup*B*dt,NHup,NH4/(B*dt))
   !if (NH4.lt.NHup*B*dt) then
   !NHup = NH4/(B*dt)
   if (NH4.lt.(NHup*B+NHr*NH4)*dt) then
      NHup = (NH4-NHr*NH4*dt)/(B*dt)
   end if
#endif

   if (QPB .le. self%maxQPB) then
      POup=self%POupmax*fT*PO4/(self%KPO4+PO4)*(1._rk-QPB/self%maxQPB)
   else
      POup=0.0_rk
   end if

#if _USE_DT_
   !POup=if(PO4.ge.POup*B*dt,POup,PO4/(B*dt))
   if (PO4.lt.POup*B*dt) then
      POup = PO4/(B*dt)
   end if
#endif

!  Avoid negative oxygene
#if _USE_DT_
   dO2= B*(self%qOB*my+self%qONO*NOup)-self%qONH*NHr*NH4-self%qOC*Cr*C
   if (O2.lt.-dO2*dt) then
!     my=0.
      NHr=0.0_rk
      Cr=0.0_rk
   end if
#endif


!  Set temporal derivatives
   _SET_ODE_(self%id_B,0*((my-self%G)*B))
   _SET_ODE_(self%id_N,0*((NOup+NHup)*B-self%G*N))
   _SET_ODE_(self%id_P,0*(POup*B-self%G*P))
   _SET_ODE_(self%id_C,0*((1-self%gam)*self%G*B-Cr*C))
   _SET_ODE_(self%id_M,0*((1-self%gam)*self%G*N-Mr*M))
   _SET_ODE_(self%id_W,0*((1-self%gam)*self%G*P-Wr*C))
   _SET_ODE_(self%id_NO3,0*(-NOup*B+NHr*NH4))
   _SET_ODE_(self%id_NH4,0*(-NHup*B-NHr*NH4+self%e*self%gam*self%G*N+Mr*M))
   _SET_ODE_(self%id_PO4,0*(-POup*B+Wr*C+self%e*self%gam*self%G*P))
   _SET_ODE_(self%id_O2,0*(B*(self%qOB*my+self%qONO*NOup)-self%qONH*NHr*NH4-self%qOC*Cr*C))

!  Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_au_daneco), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                     :: S, X, C, N
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   _GET_(self%id_salt,S)  ! local salinity

!  Retrieve current (local) state variable values.
   _GET_(self%id_C,C)
   _GET_(self%id_N,N)

    X=self%qXN*N

!  Self-shading with explicit contribution from background phytoplankton concentration.
!   _SET_EXTINCTION_(self%kw - self%eS*S + self%eC*C + self%eX*X)
   _SET_EXTINCTION_(self%kw)

!  Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   subroutine get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   class (type_au_daneco), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                     :: N, M, NO3, NH4
   real(rk)                     :: P, W, PO4
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

write(0,*) 'get_conserved_quantities'

   ! Retrieve current (local) state variable values.
   _GET_(self%id_N,N)
   _GET_(self%id_M,M)
   _GET_(self%id_NO3,NO3)
   _GET_(self%id_NH4,NH4)
   _GET_(self%id_P,P)
   _GET_(self%id_W,W)
   _GET_(self%id_PO4,PO4)

   _SET_CONSERVED_QUANTITY_(self%id_totN,N+M+NO3+NH4)
   _SET_CONSERVED_QUANTITY_(self%id_totP,P+W+PO4)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!
! !INPUT PARAMETERS:
   class (type_au_daneco),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
! !LOCAL VARIABLES:
   real(rk) :: B,N,P,C,M,W,NO3,NH4,PO4,O2 !local state variabels
   real(rk) :: NOup !Nitrate uptake
   real(rk) :: NHup !Ammonium uptake
   real(rk) :: POup !Phosphorous uptake
   real(rk) :: my   !Growth
   real(rk) :: Cr   !Carbon remineralisation
   real(rk) :: Mr   !Nitrogen remineralisation
   real(rk) :: Wr   !Phosphorous remineralisation
   real(rk) :: NHr  !Ammonium remineralisation

   real(rk) :: par,I_0,temp !local denpendency variables
   real(rk) :: my1  !Growth light
   real(rk) :: my2  !Growth nitrogen
   real(rk) :: my3  !Growth phosphorous
   real(rk) :: QNB  !Plankton nitrogen to carbon ratio
   real(rk) :: QPB  !Plankton phosphorous to carbon ratio
   real(rk) :: fT   !Temperature function
   real(rk) :: X    !Chlorophyll
   real(rk) :: QMC  !Detritus nitrogen to carbon ratio
   real(rk) :: QWC  !Detritus phosphorous to carbon ratio
!EOP
!-----------------------------------------------------------------------
!BOC
!  Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_


!  ! Assign destruction rates to different elements of the destruction matrix.
!  ! By assigning with _SET_DD_SYM_(i,j,val) as opposed to _SET_DD_(i,j,val),
!  ! assignments to dd(i,j) are automatically assigned to pp(j,i) as well.
!   _SET_DD_(self%id_B,self%id_B,self%G)
!   _SET_PP_(self%id_B,self%id_B,my)
!
!   _SET_DD_(self%id_N,self%id_N,self%G)
!   _SET_PP_(self%id_N,self%id_B,NOup+NHup)
!
!   _SET_DD_(self%id_P,self%id_P,self%G)
!   _SET_PP_(self%id_P,self%id_B,POup)
!
!   _SET_DD_(self%id_C,self%id_C,Cr)
!   _SET_PP_(self%id_C,self%id_B,(1-self%gam)*self%G)
!
!   _SET_DD_(self%id_M,self%id_M,Mr)
!   _SET_PP_(self%id_M,self%id_N,(1-self%gam)*self%G)
!
!   _SET_DD_(self%id_W,self%id_C,Wr)
!   _SET_PP_(self%id_W,self%id_P,(1-self%gam)*self%G)
!
!   _SET_DD_(self%id_NH4,self%id_B,NHup)
!   _SET_DD_(self%id_NH4,self%id_NH4,NHr)
!   _SET_PP_(self%id_NH4,self%id_N,self%e*self%gam*self%G)
!   _SET_PP_(self%id_NH4,self%id_M,Mr)
!
!   _SET_DD_(self%id_NO3,self%id_B,NOup)
!   _SET_PP_(self%id_NO3,self%id_NH4,NHr)
!
!   _SET_DD_(self%id_PO4,self%id_B,POup)
!   _SET_PP_(self%id_PO4,self%id_P,self%e*self%gam*self%G)
!   _SET_PP_(self%id_PO4,self%id_C,Wr)
!
!   _SET_DD_(self%id_O2,self%id_C,self%qOC*Cr)
!   _SET_DD_(self%id_O2,self%id_NH4,self%qONH*NHr)
!   _SET_PP_(self%id_O2,self%id_B,self%qOB*my+self%qONO*NOup)
!
!   Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine do_ppdd
!EOC

   end module fabm_au_daneco

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#endif

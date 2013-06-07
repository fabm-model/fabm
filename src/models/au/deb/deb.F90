#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_au_deb --- Dynamic Energy Budget (DEB) model for blue mussels,
! adapted for FABM by Marie Maar (MAM)
!
! !INTERFACE:
   module fabm_au_deb
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types
   use fabm_driver

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public au_deb_create
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!   Original author of au_deb: Marie Maar
!!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_au_deb
!     Variable identifiers
      type (type_state_variable_id)      :: id_volmus  !structural volume (cm3)
      type (type_state_variable_id)      :: id_rdcmus  !reserve density C  (mol-C/cm3)
      type (type_state_variable_id)      :: id_nbmus   !mussel biomass N  (mol-N/ind)
      type (type_state_variable_id)      :: id_pbmus   !mussel biomass P  (mol-P/ind)
      type (type_state_variable_id)      :: id_rpcmus  !reproductive tissue C (mol-C/ind)
      type (type_state_variable_id)      :: id_spcmus  !spawned biomass C (mol-C/m3)
      type (type_state_variable_id)      :: id_spnmus  !spawned biomass N (mol-N/m3)
      type (type_state_variable_id)      :: id_sppmus  !spawned biomass P (mol-P/m3)
      type (type_state_variable_id)      :: id_nmus    !number of mussels per m3
      type (type_state_variable_id)      :: id_decmus  !biomass of dead mussels C (mol-C/m3)
      type (type_state_variable_id)      :: id_denmus  !biomass of dead mussels N (mol-N/m3)
      type (type_state_variable_id)      :: id_depmus  !biomass of dead mussels P (mol-P/m3)
      type (type_state_variable_id)      :: id_hacmus  !biomass of harvested mussels C (mol-C/m3)
      type (type_state_variable_id)      :: id_hanmus  !biomass of harvested mussels N (mol-N/m3)
      type (type_state_variable_id)      :: id_hapmus  !biomass of harvested mussels P (mol-P/m3)
!      type (type_state_variable_id)      :: id_prey    !phytoplankton biomass from NPZD (mmol-N/m3)
      type (type_state_variable_id)      :: id_timemus    !test of time
      type (type_state_variable_id)      :: id_B       !microplankton C biomass from daneco
      type (type_state_variable_id)      :: id_C       !detritus C biomass from daneco
      type (type_state_variable_id)      :: id_N       !microplankton N biomass from daneco
      type (type_state_variable_id)      :: id_M       !detritus N biomass from daneco
      type (type_state_variable_id)      :: id_P       !microplankton P biomass from daneco
      type (type_state_variable_id)      :: id_W       !detritus P biomass from daneco
      type (type_state_variable_id)      :: id_O2      !oxygen from daneco
      type (type_state_variable_id)      :: id_NH4     !NH4 in daneco
      type (type_state_variable_id)      :: id_PO4     !PO4 in daneco
      
      type (type_diagnostic_variable_id) :: id_cowmus  !core dry weight (g/ind)
      type (type_diagnostic_variable_id) :: id_rdwmus  !reserve dry weight (g/ind)
      type (type_diagnostic_variable_id) :: id_rpwmus  !reproductive dry weight (g/ind)
      type (type_diagnostic_variable_id) :: id_bwmus   !body dry weight (g/ind)
      type (type_diagnostic_variable_id) :: id_pwmus   !population dry weight (g/ind)
      type (type_diagnostic_variable_id) :: id_lenmus  !shell length (cm)
      type (type_diagnostic_variable_id) :: id_resmus  !respiration rate (mol-O2/ind/s)
      type (type_diagnostic_variable_id) :: id_excmus  !excretion rate (mol-N/ind/s)
      type (type_dependency_id)          :: id_temp, id_salt    !celcius
      type (type_global_dependency_id)   :: id_yearday

      type (type_conserved_quantity_id)  :: id_totNmus
      type (type_conserved_quantity_id)  :: id_totPmus   

!     Model parameters
      real(rk) :: O2MIN,PIMUS,AEMUS,FH,PMAIN,VMUS,KAPPA,EG,EMAX,DVM
      real(rk) :: WMV,WME,WMR,WMF,AQDV,AQDE,AQDR,AQDF,VOLP,TL,TH,TAL,TAH,SPAWT,SPAWR
      real(rk) :: SPAWMIN,SHAPEMUS,RQC,MINC,RETZOO,OXYK
      
      contains
      
      procedure :: do
      procedure :: do_ppdd
      procedure :: get_conserved_quantities
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the deb model
!
! !INTERFACE:
   function au_deb_create(configunit,name,parent) result(self)
!
! !DESCRIPTION:
!  Here, the au namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   integer,                          intent(in)    :: configunit
   character(len=*),                 intent(in)    :: name
   class (type_model_info),target,   intent(inout) :: parent
   class (type_au_deb), pointer                    :: self
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                  :: nmus_initial= 1.0
   real(rk)                  :: volmus_initial= 0.81
   real(rk)                  :: rdcmus_initial= 0.016
   real(rk)                  :: nbmus_initial=0.002692
   real(rk)                  :: pbmus_initial=0.000168
   real(rk)                  :: rpcmus_initial=0.00
   real(rk)                  :: spcmus_initial=0.0
   real(rk)                  :: spnmus_initial=0.0
   real(rk)                  :: sppmus_initial=0.0
   real(rk)                  :: decmus_initial=0.0
   real(rk)                  :: denmus_initial=0.0
   real(rk)                  :: depmus_initial=0.0
   real(rk)                  :: hacmus_initial=0.0
   real(rk)                  :: hanmus_initial=0.0
   real(rk)                  :: hapmus_initial=0.0
   real(rk)                  :: timemus_initial=1.0
 
   real(rk)                  :: O2MIN = 90.          !Oxygen threshold (mmol/m3)
   real(rk)                  :: PIMUS = 0.0038       !max. ingestion (mol-C/cm2/)
   real(rk)                  :: AEMUS = 0.75         !assimilation efficiency
   real(rk)                  :: FH = 0.016           !half saturation constant (mol-C/m3)
   real(rk)                  :: PMAIN = 4.630E-04    !Maintenance (mol-C/cm3/d)
   real(rk)                  :: VMUS = 80.5          !max.volume (cm3) 4.31 cm
   real(rk)                  :: KAPPA = 0.7          !fraction spent on growth+maintenance
   real(rk)                  :: EG = 0.02885         !cost of growth (mol-C/cm3)
   real(rk)                  :: EMAX = 0.03326       !maximum reserve density (mol-C/cm3)
   real(rk)                  :: DVM = 2.05E-02       !density of structural tissue (mol-C/cm3)
   real(rk)                  :: WMV = 23.9           !molar weight structural tissue (g/mol)
   real(rk)                  :: WME = 28.8           !molar weight reserves (g/mol)
   real(rk)                  :: WMR = 28.8           !molar weight reproductive tissue (g/mol)
   real(rk)                  :: WMF = 24.6           !molar weight food (g/mol)
   real(rk)                  :: AQDV = 0.70          !aqueous fraction
   real(rk)                  :: AQDE = 0.80          !aqueous fraction
   real(rk)                  :: AQDR = 0.80          !aqueous fraction
   real(rk)                  :: AQDF = 0.75          !aqueous fraction
   real(rk)                  :: VOLP = 0.06          !volume at puperty(cm3)
   real(rk)                  :: TL = 275.            !lower temp. boundary (K)
   real(rk)                  :: TH = 296.            !upper temp. boundary (K)
   real(rk)                  :: TAL = 45430.         !Arrhenius temp. lower (K)
   real(rk)                  :: TAH = 31376.         !Arrhenius temp. upper (K)
   real(rk)                  :: SPAWT = 15.          !min. spawning temp. (Celcius)
   real(rk)                  :: SPAWR = 0.523        !spawning rate (/d)
   real(rk)                  :: SPAWMIN = 0.30       !minimum gonadosomaticratio
   real(rk)                  :: SHAPEMUS = 0.314     !0.287 shape factor
   real(rk)                  :: RQC = 0.84           !respiratory qoutient
   real(rk)                  :: MINC = 3.            !filtration threshold (mmol/m3)
   real(rk)                  :: RETZOO = 0.3         !retention eff. of zooplankton
   real(rk)                  :: OXYK= 25.            !half sat. O2 (mmol/m3)
!   character(len=64)         :: prey_source_variable =''  !input var. from NPZD
   character(len=64)         :: B_source_variable =''  !input var. from Daneco
   character(len=64)         :: C_source_variable =''  !input var. from Daneco
   character(len=64)         :: P_source_variable =''  !input var. from Daneco
   character(len=64)         :: N_source_variable =''  !input var. from Daneco
   character(len=64)         :: O2_source_variable =''  !input var. from Daneco
   character(len=64)         :: M_source_variable =''  !input var. from Daneco
   character(len=64)         :: W_source_variable =''  !input var. from Daneco
   character(len=64)         :: NH4_target_variable =''  !output var. to Daneco
   character(len=64)         :: PO4_target_variable =''  !output var. to Daneco
   
   real(rk), parameter :: secs_pr_day = 86400.
   namelist /au_deb/ volmus_initial,rdcmus_initial, nbmus_initial, rpcmus_initial, &
                  spcmus_initial, spnmus_initial, sppmus_initial, nmus_initial,    &
                  decmus_initial, denmus_initial, depmus_initial, hacmus_initial,  &
                  hanmus_initial,pbmus_initial, hapmus_initial, timemus_initial,   &
                  O2MIN,PIMUS,AEMUS,FH,PMAIN,VMUS,                                 &
                  KAPPA,EG,EMAX,DVM,WMV,WME,WMR,WMF,AQDV,AQDE,AQDR,                &
                  AQDF,VOLP,TL,TH,TAL,TAH,SPAWT,SPAWR,SPAWMIN,SHAPEMUS,RQC,        &
                  MINC,RETZOO,OXYK,                                                &
                  B_source_variable, C_source_variable, O2_source_variable,        &
                  N_source_variable, M_source_variable,                            &
                  NH4_target_variable,PO4_target_variable,                         &
                  P_source_variable, W_source_variable
!EOP 
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(configunit,nml=au_deb,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%O2MIN      = O2MIN
   self%PIMUS      = PIMUS/secs_pr_day
   self%AEMUS      = AEMUS
   self%FH         = FH
   self%PMAIN      = PMAIN/secs_pr_day
   self%VMUS       = VMUS
   self%KAPPA      = KAPPA
   self%EG         = EG
   self%EMAX       = EMAX
   self%DVM        = DVM
   self%WMV        = WMV
   self%WME        = WME
   self%WMR        = WMR
   self%WMF        = WMF
   self%AQDV       = AQDV
   self%AQDE       = AQDE
   self%AQDR       = AQDR
   self%AQDF       = AQDF
   self%VOLP       = VOLP
   self%TL         = TL
   self%TH         = TH
   self%TAL        = TAL
   self%TAH        = TAH
   self%SPAWT      = SPAWT
   self%SPAWR      = SPAWR/secs_pr_day
   self%SPAWMIN    = SPAWMIN
   self%SHAPEMUS   = SHAPEMUS
   self%RQC        = RQC
   self%MINC       = MINC
   self%RETZOO     = RETZOO
   self%OXYK       = OXYK

   
!  Register state variables
!  Pelagic
   call self%register_state_variable(self%id_volmus,'volmus','cm**3','structural volume',         &
                                     volmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_rdcmus,'rdcmus','mol-C/ind','reserve density C',     &
                                     rdcmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_nbmus,'nbmus','mol-N/ind','mussel biomass N',        &
                                     nbmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_pbmus,'pbmus','mol-P/ind','mussel biomass P',        &
                                     pbmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_rpcmus,'rpcmus','mol-C/ind','reproductive tissue C', &
                                     rpcmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_spcmus,'spcmus','mol-C/ind','spawned biomass C',     &
                                     spcmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_spnmus,'spnmus','mol-N/ind','spawned biomass N',     &
                                     spnmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_sppmus,'sppmus','mol-P/ind','spawned biomass P',     &
                                     sppmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_nmus,'nmus','ind/m3','number of individuals',      &
                                     nmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_decmus,'decmus','mol-C/m2','dead biomass C',         &
                                     decmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_denmus,'denmus','mol-N/m2','dead biomass N',         &
                                     denmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_depmus,'depmus','mol-P/m2','dead biomass P',         &
                                     depmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_hacmus,'hacmus','mol-C/m2','harvest biomass C',      &
                                     hacmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_hanmus,'hanmus','mol-N/m2','harvest biomass N',      &
                                     hanmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_hapmus,'hapmus','mol-P/m2','harvest biomass P',      &
                                     hapmus_initial,minimum=_ZERO_,no_river_dilution=.true.)
   call self%register_state_variable(self%id_timemus,'timemus','day','time',      &
                                     timemus_initial,minimum=_ZERO_,no_river_dilution=.true.)

   ! Register link to external pelagic prey and mineral pools.
   ! Prey will be used to feed upon, mineral pool to place waste products in.
!   self%id_prey = self%register_state_dependency(prey_source_variable)
   call self%register_state_dependency(self%id_B,  B_source_variable)
   call self%register_state_dependency(self%id_C,  C_source_variable)
   call self%register_state_dependency(self%id_P,  P_source_variable)
   call self%register_state_dependency(self%id_N,  N_source_variable)
   call self%register_state_dependency(self%id_O2, O2_source_variable)
   call self%register_state_dependency(self%id_M,  M_source_variable)
   call self%register_state_dependency(self%id_W,  W_source_variable)
   call self%register_state_dependency(self%id_NH4,NH4_target_variable)
   call self%register_state_dependency(self%id_PO4,PO4_target_variable)

!  Register diagnostic variables

   call self%register_diagnostic_variable(self%id_cowmus,'cowmus','g-DW/ind','core dry weigh',          &
                                          time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_rdwmus,'rdwmus','g-DW/ind', 'reserve dry weight',     &
                                          time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_rpwmus,'rpwmus','gDW/ind', 'reproductive dry weight', &
                                          time_treatment=time_treatment_step_integrated)

   call self%register_diagnostic_variable(self%id_bwmus,'bwmus','g-DW/ind','dry weight mussel',         &
                                          time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_pwmus,'pwmus','g-DW/m3','dry weight cohort',          &
                                          time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_lenmus,'lenmus','cm','shell length',                  &
                                          time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_resmus,'resmus','mol-O2/ind/s','respiration',         &
                                          time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_excmus,'excmus','mol-N/ind/s','excretion',            &
                                          time_treatment=time_treatment_step_integrated)


  ! Register environmental dependencies
   call self%register_dependency(self%id_temp,   varname_temp)
   call self%register_dependency(self%id_salt,   varname_salt)
   call self%register_dependency(self%id_yearday,varname_yearday)

   ! Register conserved quantities
   call self%register_conserved_quantity(self%id_totNmus,'N','mmol/m**3','totNmus')
   call self%register_conserved_quantity(self%id_totPmus,'P','mmol/m**3','totPmus')
   
   return

99 call fatal_error('au_deb_create','Error reading namelist au_deb')
   
   end function au_deb_create
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_au_deb), intent(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding

! !LOCAL VARIABLES:
   real(rk)                   :: nmus, volmus,rdcmus, nbmus, pbmus
   real(rk)                   :: rpcmus, spcmus, spnmus, sppmus
   real(rk)                   :: decmus, denmus, depmus, hacmus, hanmus, hapmus
   real(rk)                   :: lenmus,cowmus,rdwmus,bwmus,pwmus
   real(rk)                   :: rpwmus,cocmus,recmus,cbmus, dwmus
   real(rk)                   :: temp, salt, yearday
!   real(rk)                   :: prey
   real(rk), parameter        :: secs_pr_day = 86400.
   real(rk)                   :: TMP
   real(rk)                   :: XING
   real(rk)                   :: FOOD
   real(rk)                   :: AFC, PAMUS
   real(rk)                   :: FACMUS, FANMUS, FAPMUS
   real(rk)                   :: PHYI, ZOOI, DETI
   real(rk)                   :: MRESP
   real(rk)                   :: RESP, SPAW
   real(rk)                   :: RDEN, QRES, RDENN, RDENP
   real(rk)                   :: CATAB, STARV
   real(rk)                   :: VOLG, EXMUS, LMUS
   real(rk)                   :: REPRO, OXYM, HARV
   real(rk)                   :: QPRES, resmus, excmus, dt
   real(rk)                   :: B, C, N, M, O2, NH4, P, W, PO4
   real(rk)                   :: MPC, ZPC, DET, DAY, timemus, temp1
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_
!KB
!  Retrieve current (local) state variable values.
   _GET_(self%id_nmus,nmus)
   _GET_(self%id_volmus,volmus)
   _GET_(self%id_rdcmus,rdcmus)
   _GET_(self%id_nbmus,nbmus)
   _GET_(self%id_pbmus,pbmus)
   _GET_(self%id_rpcmus,rpcmus)
   _GET_(self%id_spcmus,spcmus)
   _GET_(self%id_spnmus,spnmus)
   _GET_(self%id_spnmus,sppmus)
   _GET_(self%id_decmus,decmus)
   _GET_(self%id_denmus,denmus)
   _GET_(self%id_depmus,depmus)
   _GET_(self%id_hacmus,hacmus)
   _GET_(self%id_hanmus,hanmus)
   _GET_(self%id_hapmus,hapmus)
!   _GET_(self%id_prey,prey)     ! prey density - pelagic NPZD
   _GET_(self%id_timemus,timemus)     
   _GET_(self%id_B,B)   
   _GET_(self%id_C,C)   
   _GET_(self%id_N,N)   
   _GET_(self%id_O2,O2)  
   _GET_(self%id_M,M)   
   _GET_(self%id_P,P)
   _GET_(self%id_W,W)     

 ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)  ! local water temperature
   _GET_(self%id_salt,salt)  ! local salinity
   _GET_GLOBAL_(self%id_yearday,yearday)

 ! Define some intermediate quantities that will be reused multiple times.

 !input variables from data for testing
!   temp1 = 8.7670E-04*timemus**2. - 6.8592E-01*timemus +135.5285
!   MPC = (2.438+26.63*exp(-0.5*((timemus-434.82)/9.86)**2))*50./12.*1.1 !prey*6.625*20. !microzooplankton biomass
!   ZPC = 0.d0 !zoo
!   DET = 0.d0  !detritus
!   TOTFOOD = MPC + ZPC*self%RETZOO + DET
!   O2 = 300. !oxygen

 !input variables from Daneco
   dt = 60.0
   MPC = (B + C)
   ZPC = 0.d0
! 
!---core, reserv and total carbon biomass [mol-C/ind]
   cocmus = volmus*self%DVM
   recmus = rdcmus*volmus
   cbmus = cocmus + recmus + rpcmus
!---mussel dry weight biomass [g-DW/ind] (Ren & Ross 2005)
   dwmus = cocmus*self%WMV*(1.0-self%AQDV) + &
           recmus*self%WME*(1.0-self%AQDE) + rpcmus*self%WMR*(1.0-self%AQDR)
!---shell length
   LMUS = volmus**0.333/self%SHAPEMUS
!---temperature function mussel [0:1] 
   TMP = TEMMUS(self,temp1)
!---mussel food saturation function
   FOOD = FOODMUS(self,MPC,ZPC) 
!---mussel total ingestion [mol-C/ind/s]
   XING = INGMUS(self,TMP,FOOD,volmus,MPC,O2)
!---ingestion microplankton from Daneco
  if ((B*0.001.gt.XING*B/MPC*nmus*dt).and.(N*0.001.gt.XING*N/MPC*nmus*dt).and.(P*0.001.gt.XING*P/MPC*nmus*dt)) then
   PHYI = XING*B/MPC
  else
   PHYI = 0.0
  end if
!---ingestion detritus from daneco
  if ((C*0.001.gt.XING*C/MPC*nmus*dt).and.(M*0.001.gt.XING*M/MPC*nmus*dt).and.(W*0.001.gt.XING*W/MPC*nmus*dt)) then
   DETI = XING*C/MPC
  else
   DETI = 0.0
  end if
!---corrected ingestion (mol-C/ind/s)
   XING = PHYI + DETI
!---maximum assimilation (mol-C/ind/s)
   PAMUS = self%PIMUS*self%AEMUS
!---mussel assimilation [mol-C/ind/s]
   AFC = XING*self%AEMUS
!---maintenance [mol-C/ind/s]
   MRESP = MAINMUS(self,TMP,volmus)
!---energy density [mol-C/cm3/s]
   RDEN = DENSMUS(self,FOOD,rdcmus,volmus,MPC,PAMUS,TMP)
!---catabolic rate [mol-C/ind/s]
   CATAB = CATABMUS(self,MRESP,rdcmus,volmus,PAMUS,TMP)
!---volume growth [cm3/s]
   VOLG = VOLMMUS(self,MRESP,CATAB)
!---reproductive tissue [mol-C/ind/s]
   REPRO = REPRMUS(self,TMP,CATAB,volmus)
!---spawned biomass [mol-C/ind/s]
   SPAW = SPAWNING(self,rpcmus,cocmus,temp1)
!---reprocductive tissue pay for maintenance if VOLG is neg.
   if (VOLG.lt.0.0) then
     REPRO = REPRO + VOLG*self%DVM
     VOLG  = 0.0
   end if
   
  if (rpcmus.lt.0.0) then
     denmus = nmus*nbmus
     decmus = nmus*cbmus
     depmus = nmus*pbmus
     REPRO = 0.0
     nmus= 0.0
     volmus = 0.0
     rpcmus = 0.0
     rdcmus= 0.0
     nbmus = 0.0
     pbmus = 0.0
     cbmus = 0.0 
  endif

!  Coupling (output) to Daneco model /remember conversion to mmol
!---mussel defecation C [mol-C/m3/s]
   FACMUS = XING-AFC
!---mussel defecation N [mol-N/m3/s]
   FANMUS = (PHYI*N/B + DETI*M/C)*(1.0 - self%AEMUS)
!---mussel defecation N [mol-P/m3/s]
   FAPMUS = (PHYI*P/B + DETI*W/C)*(1.0 - self%AEMUS)
!---respiration [mol-02/ind/s]
   RESP = RESPMUS(self,CATAB,VOLG,REPRO)
!---ingestion ZPC
!   ZOOI = XING*ZPC*self%RETZOO/TOTFOOD
!---starvartion mortality
   STARV = MORTALITY(self,rdcmus,volmus,nmus)
!---oxygen mortality
   OXYM = OXYMORT(self,O2,nmus)
!---harvesting
   if (yearday.eq.90.0) then
     HARV = HARVESTING(self,LMUS,nmus)
   end if
!---N:C ratio mussel biomass
   QRES = nbmus/cbmus
!---P:N ratio mussel biomass
   QPRES = pbmus/nbmus
!---mussel excretion of N [mol-N/ind/s]
   EXMUS = RESP*QRES/self%RQC
!---nitrogen growth [mol-N/ind/s]
   RDENN = (PHYI*N/B + DETI*M/C)*self%AEMUS - EXMUS
!---phosphorous growth [mol-P/ind/s]
   RDENP = (PHYI*P/B + DETI*W/C)*self%AEMUS - EXMUS*QPRES
   
   ! Set temporal derivatives
   _SET_ODE_(self%id_volmus, VOLG)
   _SET_ODE_(self%id_nmus,(-STARV - OXYM -HARV)*nmus)
   _SET_ODE_(self%id_rdcmus,RDEN)
   _SET_ODE_(self%id_nbmus,RDENN)
   _SET_ODE_(self%id_pbmus,RDENP)
   _SET_ODE_(self%id_rpcmus,REPRO-SPAW)
   _SET_ODE_(self%id_spcmus,SPAW*nmus*cbmus)
   _SET_ODE_(self%id_spnmus,SPAW*nmus*nbmus)
   _SET_ODE_(self%id_sppmus,SPAW*nmus*pbmus)
   _SET_ODE_(self%id_decmus,(STARV+OXYM)*nmus*cbmus)
   _SET_ODE_(self%id_denmus,(STARV+OXYM)*nmus*nbmus)
   _SET_ODE_(self%id_depmus,(STARV+OXYM)*nmus*pbmus)
   _SET_ODE_(self%id_hacmus,HARV*nmus*cbmus)
   _SET_ODE_(self%id_hanmus,HARV*nmus*nbmus)
   _SET_ODE_(self%id_hapmus,HARV*nmus*pbmus)
   _SET_ODE_(self%id_timemus,1./86400.)

   ! link to daneco variables
   _SET_ODE_(self%id_B, -PHYI*nmus*1000.d0)
   _SET_ODE_(self%id_N, -PHYI*N/B*nmus*1000.d0)
   _SET_ODE_(self%id_P, -PHYI*P/B*nmus*1000.d0)
   _SET_ODE_(self%id_C, (-DETI + FACMUS)*nmus*1000.d0)
   _SET_ODE_(self%id_M, (-DETI*M/C + FANMUS)*nmus*1000.d0)
   _SET_ODE_(self%id_W, (-DETI*W/C + FAPMUS)*nmus*1000.d0)
   _SET_ODE_(self%id_O2, -RESP*nmus*1000.d0)
   _SET_ODE_(self%id_NH4, EXMUS*nmus*1000.d0)
   _SET_ODE_(self%id_PO4, EXMUS*QPRES*nmus*1000.d0)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_cowmus,volmus*self%WMV*self%DVM*(1.0-self%AQDV))
   _SET_DIAGNOSTIC_(self%id_rdwmus,rdcmus*volmus*self%WME*(1.0-self%AQDE))
   _SET_DIAGNOSTIC_(self%id_rpwmus,rpcmus*self%WMR*(1.0-self%AQDR))
   _SET_DIAGNOSTIC_(self%id_pwmus,dwmus*nmus)
   _SET_DIAGNOSTIC_(self%id_bwmus,dwmus) 
   _SET_DIAGNOSTIC_(self%id_lenmus,LMUS)
   _SET_DIAGNOSTIC_(self%id_resmus,RESP)
   _SET_DIAGNOSTIC_(self%id_excmus,EXMUS)
   
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine do
!EOC
!----------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   class (type_au_deb),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:

   end subroutine do_ppdd
!EOC
!
!--------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ARRHENIUS TEMPERATURE FUNCTION FOR MUSSELS  [-]
!      
!  !INTERFACE:
    pure real(rk) function TEMMUS(self,temp1)
!   
! !DESCRIPTION:
! 
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
!
   type (type_au_deb), intent(in)   :: self
   real(rk), intent(in)             :: temp1
   real(rk)                         :: TMP1
! 
! !REVISION HISTORY:
!  Original author(s): MAM   
!  after van der Meer (2006), J Sea Res 56:107-224 
!  Kooijman 2000 Dynamic energy and mass budgets 
!
!EOP
!----------------------------------------------------------------
!BOP
   TMP1 = temp1 + 273.15

   TEMMUS=(1.+EXP(self%TAL/TMP1-self%TAL/self%TL)+EXP(self%TAH/self%TH-self%TAH/TMP1))**-1.0

   RETURN

   end function TEMMUS
!EOP
!---------------------------------------------------------------------------
!BOP
! !IROUTINE: FOOD SATURATION FUNCTION FOR MUSSEL INGESTION [-]
!      
!  !INTERFACE:
   pure real(rk) function FOODMUS(self,MPC,ZPC) 
! 
! !DESCRIPTION:
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in)   :: self
  real(rk) , intent(in)            :: MPC,ZPC
  real(rk)                         :: PHY

! !REVISION HISTORY:
!  Original author(s): MAM   
!  After van der Meer 2006, p.86 and table 3
!EOP
!-------------------------------------------------------------------
!BOP
           
     PHY = (MPC + ZPC*self%RETZOO)*0.001

     FOODMUS = PHY/(PHY+self%FH)

     RETURN

     end function FOODMUS
!EOP
!------------------------------------------------------------------------
!BOP
! !IROUTINE: Ingestion rate of mussels [mol-C/ind/s]
!      
!  !INTERFACE:
   pure real(rk) function INGMUS(self,TMP,FOOD,volmus, MPC, O2)
! 
! !DESCRIPTION:
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: TMP, volmus, FOOD, MPC, O2
  real(rk)                       :: OXY

! !REVISION HISTORY:
!  Original author(s): MAM   
!  After van der Meer 2006, p.86 
!EOP    
!-------------------------------------------------------------------
!BOP

     IF (O2.LT.self%O2MIN) THEN
       OXY= O2**3/(O2**3+ self%OXYK**3)
     ELSE
       OXY= 1.0
     END if

     if (volmus.lt.self%VMUS) then
      IF (MPC.GT.self%MINC) THEN
        INGMUS =  self%PIMUS*volmus**0.666*FOOD*TMP*OXY
      ELSE
        INGMUS = 0.0
      end if
     ELSE
        INGMUS = 0.0
     end if    

     RETURN

     end function INGMUS    
!EOP
!------------------------------------------------------------------------
!BOP
! !IROUTINE: MAINTENANCE BY MUSSELS* [molC/s]
!      
!  !INTERFACE:
   pure real(rk) function MAINMUS(self,TMP,volmus)
! 
! !DESCRIPTION:
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: TMP, volmus
!
! !REVISION HISTORY:
!  Original author(s): MAM   
!  eq.4, 2.term van der Meer 2006
!EOP    
!-------------------------------------------------------------------
!BOP
     MAINMUS = self%PMAIN*volmus*TMP

     RETURN

     END function MAINMUS

!------------------------------------------------------------------------
!BOP
! !IROUTINE: 
!      
!  !INTERFACE:
   pure real(rk) function CATABMUS(self,MRESP,rdcmus,volmus,PAMUS,TMP)
! 
! !DESCRIPTION: CATABOLIC RATE [mol-C/s]
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: rdcmus, volmus, MRESP,PAMUS,TMP
  real(rk)                       :: RHS1, RHS2
!
! !REVISION HISTORY:
!  Original author(s): MAM   
!  eq.7 van der Meer 2006 
!EOP    
!-------------------------------------------------------------------
!BOP
     RHS1 = rdcmus/(self%KAPPA*rdcmus+self%EG)
     RHS2 = PAMUS*TMP*self%EG/self%EMAX*volmus**0.666 + MRESP
       
     CATABMUS = RHS1*RHS2
           
     RETURN

     END function CATABMUS
!------------------------------------------------------------------------
!BOP
! !IROUTINE: VOLUME GROWTH [cm3/s]
!      
!  !INTERFACE:
   pure real(rk) function VOLMMUS(self,MRESP,CATAB)
! 
! !DESCRIPTION:
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: MRESP, CATAB
!
! !REVISION HISTORY:
!  Original author(s): MAM   
!  eq.4, van der Meer 2006
!  DW:p.352. Ren and Ross 2005/table 2 van der Veer 
!EOP    
!-------------------------------------------------------------------
!BOP
     VOLMMUS = (CATAB*self%KAPPA-MRESP)/self%EG
        
     RETURN

     END function VOLMMUS
     
!------------------------------------------------------------------------
!BOP
! !IROUTINE: ENERGY DENSITY [molC/cm3/s] 
!      
!  !INTERFACE:
   pure real(rk) function DENSMUS(self,FOOD,rdcmus,volmus,MPC,PAMUS,TMP)
! 
! !DESCRIPTION:
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: FOOD,rdcmus,volmus,MPC,PAMUS,TMP
  real(rk)                       :: RHS1
!
! !REVISION HISTORY:
!  Original author(s): MAM   
!  eq.2 van der Meer 2006
!EOP    
!-------------------------------------------------------------------
!BOP

     RHS1 = (PAMUS/volmus**0.333)*TMP
    
     IF (MPC.GE.self%MINC) THEN
        DENSMUS = RHS1*(FOOD-rdcmus/self%EMAX)
     ELSE
        DENSMUS = RHS1*(-rdcmus/self%EMAX)         
     END IF
     
     RETURN

     END function DENSMUS
     
!------------------------------------------------------------------------
!BOP
! !IROUTINE: 
!      
!  !INTERFACE:
   pure real(rk) function RESPMUS(self,CATAB, VOLG, REPRO)
! 
! !DESCRIPTION:
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: CATAB, REPRO, VOLG
!
! !REVISION HISTORY:
!  Original author(s): MAM      
!EOP    
!-------------------------------------------------------------------
!BOP

     RESPMUS = self%RQC*(CATAB - REPRO - VOLG*self%DVM)


     RETURN

     END function RESPMUS
!EOP    
!-------------------------------------------------------------------
!BOP
! !IROUTINE: REPRODUCTIVE TISSUE [mol-C/ind/s]
!      
!  !INTERFACE:
   pure real(rk) function REPRMUS(self,TMP,CATAB,volmus)
! 
! !DESCRIPTION:
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: TMP, CATAB, volmus
  real(rk)                       :: RHS1, PJMUS
!
! !REVISION HISTORY:
!  Original author(s): MAM   
!  eq.9 van der Meer 2006   
!EOP    
!-------------------------------------------------------------------
!BOP
     RHS1 = MIN(volmus,self%VOLP)*self%PMAIN*TMP
     PJMUS = RHS1*(1.-self%KAPPA)/self%KAPPA
     REPRMUS = (1.-self%KAPPA)*CATAB-PJMUS      
         
     RETURN

     END function REPRMUS
     
!------------------------------------------------------------------------
!BOP
! !IROUTINE: SPAWNING (mol-C/ind/s)
!      
!  !INTERFACE:
   pure real(rk) function SPAWNING(self,rpcmus,cocmus,temp1)
! !DESCRIPTION:
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: temp1,rpcmus,cocmus
  real(rk)                       :: GONADO
!-------------------------------------------------------------------
!BOP
     GONADO = rpcmus/(rpcmus + cocmus)
     IF ((GONADO.GT.self%SPAWMIN).AND.(temp1.GT.self%SPAWT)) THEN
          SPAWNING =self%SPAWR*rpcmus
     else
          SPAWNING = 0.0
                   
     END IF

     RETURN

     end function SPAWNING
!EOP
!------------------------------------------------------------------------

!BOP
! !IROUTINE: Starvation mortality  [/s]
!      
!  !INTERFACE:
   pure real(rk) function MORTALITY(self,rdcmus,volmus,nmus)
! 
! !DESCRIPTION:
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk),  intent(in)          :: rdcmus,volmus,nmus
  real(rk)                       :: ERATIO
!
! !REVISION HISTORY:
!  Original author(s): MAM      
!EOP    
!-------------------------------------------------------------------
!BOP

     real(rk), parameter        :: ALPHA = 70.
     real(rk), parameter        :: BETA = 100.

     ERATIO = rdcmus/self%EMAX

     if(nmus.gt.1.0) then
       MORTALITY = (1.-1./(1.+BETA*EXP(-ERATIO*ALPHA)))/86400.
     else
       MORTALITY = 0.0
     end if
 
     RETURN

     END function MORTALITY
     
    
!------------------------------------------------------------------------
!BOP
! !IROUTINE: HARVESTING OF MUSSELS
!      
!  !INTERFACE:
   pure real(rk) function HARVESTING(self,lenmus,nmus) 
! 
! !DESCRIPTION: 
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: lenmus,nmus
!
! !REVISION HISTORY:
!  Original author(s): MAM      
!EOP    
!-------------------------------------------------------------------
!BOP    
   !  if (yearday.eq.90) then
      if ((lenmus.gt.5.0).and.(nmus.gt.1.0)) then

       HARVESTING = 1.0/86400.
      else
       HARVESTING = 0.0
      end if
   !  end if
    
     RETURN
       
     END function HARVESTING
     
!------------------------------------------------------------------------
!BOP
! !IROUTINE: OXYGEN mortality  [/s]
!      
!  !INTERFACE:
   pure real(rk) function OXYMORT(self,O2,nmus)
! 
! !DESCRIPTION: 
!  
! !USES:
   implicit none
!
! INPUT PARAMETERS:   
  type (type_au_deb), intent(in) :: self
  real(rk), intent(in)           :: O2,nmus
  real(rk)                       :: OXY
!
! !REVISION HISTORY:
!  Original author(s): MAM      
!EOP    
!-------------------------------------------------------------------
!BOP  
     real(rk), parameter        :: MMAX= 0.012 !per day

     IF ((O2.LT.self%O2MIN).and.(nmus.gt.1.0)) THEN       
         OXY= O2**3*(O2**3+self%OXYK**3)**-1.
         OXYMORT = MMAX*(1.0-OXY)*(86400.)**-1.
      ELSE
        OXYMORT = 0.0
      END IF

      RETURN

      END function OXYMORT
    
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities 
!
! !INTERFACE:
   subroutine get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   class (type_au_deb), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman, Marie Maar
!
! !LOCAL VARIABLES:
   real(rk)                     :: spnmus,denmus,hanmus,nbmus,sppmus,depmus,hapmus,pbmus
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   _GET_(self%id_spnmus,spnmus)
   _GET_(self%id_denmus,denmus)
   _GET_(self%id_hanmus,hanmus)
   _GET_(self%id_nbmus, nbmus)
   _GET_(self%id_sppmus,sppmus)
   _GET_(self%id_depmus,depmus)
   _GET_(self%id_hapmus,hapmus)
   _GET_(self%id_pbmus, pbmus)

   _SET_CONSERVED_QUANTITY_(self%id_totNmus,(spnmus+denmus+hanmus+nbmus)*1000.)
   _SET_CONSERVED_QUANTITY_(self%id_totPmus,(sppmus+depmus+hapmus+pbmus)*1000.)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine get_conserved_quantities
!EOC

!-----------------------------------------------------------------------

   end module fabm_au_deb

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

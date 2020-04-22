#include "fabm_driver.h"

! -------------------------------------------------------------------------------------------------------------!
! This is a FABM (Fortran) code to model protist trophic strategies                                            !
! This model is a development of previous works (Flynn & Mitra 2009 J Plankton Res, Mitra et al. 2016 Protist) ! 
! and is described in full in Leles et al 2018 J Plankton Res.                                                 !
!                                                                                                              !
! The code requires external dependencies from ERSEM, which can be obtained upon registration at:              !
! http://ersem.com                                                                                             !
! -------------------------------------------------------------------------------------------------------------!

module su_mixo

  use fabm_types
  use fabm_particle
  use fabm_expressions

  implicit none

  private

   ! Definitions below copied from ~/ersem/src/shared.F90
#ifdef IRON
   logical,parameter :: use_iron = .true.
#else
   logical,parameter :: use_iron = .false.
#endif

   type (type_bulk_standard_variable),parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_calcite_in_biota = type_bulk_standard_variable(name='total_calcite_in_biota',units='mg C/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: photosynthesis_rate = type_bulk_standard_variable(name='photosynthesis_rate',units='mg C/m^3/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: particulate_organic_absorption_coefficient = type_bulk_standard_variable(name='particulate_organic_absorption_coefficient',units='1/m',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: particulate_organic_backscatter_coefficient = type_bulk_standard_variable(name='particulate_organic_backscatter_coefficient',units='1/m',aggregate_variable=.true.)

  type,extends(type_particle_model),public :: type_su_mixo

    ! Variables (protist and external pools)
    type(type_state_variable_id) :: id_mF, id_mFN, id_mFP, id_mFChl, id_mC, id_mN, id_mP, id_mChl                 
    type(type_state_variable_id) :: id_NH4, id_NO3, id_DIP, id_DIC, id_sDOMC, id_sDOMN, id_sDOMP, id_DOC, id_RPs, id_POC, id_PON, id_POP
    type(type_state_variable_id) :: id_mSi, id_N5s                       ! allowing the description of diatoms
    type(type_bottom_state_variable_id) :: id_bc, id_bp, id_bn, id_bsi   ! allowing sinking  
    ! Variables related to the predation function
    type (type_model_id),               allocatable,dimension(:) :: id_prey
    type (type_dependency_id),          allocatable,dimension(:) :: id_preyc,id_preyn,id_preyp,id_preys,id_preyf,id_preyl,id_preychl
    type (type_state_variable_id),      allocatable,dimension(:) :: id_preyf_target
    type (type_diagnostic_variable_id), allocatable,dimension(:) :: id_sI, id_Cp, id_I, id_Ccon, id_Icap
    ! Other diagnostics
    type(type_diagnostic_variable_id) :: id_mNC, id_mPC, id_mChlC, id_mFC, id_mFChlC, id_mFNC, id_mFPC, id_mNCu, id_mPCu, id_mNPCu, id_mupN, id_mIPu,id_PSrel, id_mPS, id_mCas, id_Rphot, id_Rhet, id_mCu, id_mDgC, id_FPS, id_mPStot, id_kpmax, id_FCmax, id_Cfix, id_Cresp, id_mSC, id_SCu, id_Sup, id_photCup, id_hetCup, id_photNup, id_hetNup, id_photPup, id_hetPup, id_PAR_PB, id_mPqm, id_Scon, id_Pbal, id_BalRes, id_PbalCon, id_PD, id_Dmax, id_FCrelA, id_AE_op, id_Urel, id_mRTot, id_et, id_mregN, id_mregP, id_Preas, id_NH4reas, id_npp, id_FCrelV, id_stress, id_mortmC, id_UmpN, id_mAVP, id_mNVP, id_mAV, id_mNV, id_mfrat, id_mfrNC, id_mupN_mC, id_mNH4up, id_mNO3up, id_mfrPC, id_mIPu_mC, id_ImF, id_ImFChl,  id_ImFk, id_mDgChl_mC, id_Chl_decay, id_Chl_not_used, id_mDOCout, id_mVOCL, id_mVONL, id_mVOPL, id_mortC, id_mortN, id_mortP 
    ! External dependencies
    type(type_dependency_id)          :: id_mCudep, id_avgCu, id_PAR, id_ETW, id_mPStotdep, id_avgPStot

    ! Parameters
    real(rk) :: AAsyn, redco
    real(rk) :: AEmax, AEmin
    real(rk) :: mBR, mMR
    real(rk) :: mDOCpc
    real(rk) :: mChlCabs, mM, Pbalcrit
    real(rk) :: FCabs, FCmin
    real(rk) :: Has, Kas
    real(rk) :: Heq, Keq
    real(rk) :: Hhet, Khet
    real(rk) :: Hing, King
    real(rk) :: Hpbal, KPbal
    real(rk) :: Hpd, Kpd
    real(rk) :: Hrc, Krc
    real(rk) :: mQh, mKxi, mbeta
    real(rk) :: mAKt, mNKt, mPKu
    real(rk) :: Kec
    real(rk) :: mKQN, KQP
    real(rk) :: NCabs, NCm, mNCm, mNCo
    real(rk) :: mPCm, PCabs, PCm, PCo
    real(rk) :: mApref, mNpref, mPpref
    integer  :: Smix, Spd, Stype, Svol
    real(rk) :: malpha, aalpha
    real(rk) :: mUm, mUphot
    real(rk) :: DChl
    real(rk) :: MetMult
    integer  :: nprey                                           
    real(rk),allocatable :: preyESD(:), pref(:)    
    real(rk) :: mESD, pi, wSCEB, aSCEB, bSCEB    
    real(rk) :: Smin, Smax, Sopt                                       
    integer ,allocatable :: steal(:)            
    real(rk) :: kpmax                                     
    real(rk) :: iopABS, iopBBS                  
    logical  :: SSi                             
    real(rk) :: KSi, mSCm, SCo, betaSi, SCabs   
    real(rk) :: q10                             
    real(rk) :: sd                              
    real(rk) :: diss_frac                        
    real(rk) :: sink
  contains
    ! Add model procedures here
    procedure :: initialize
    procedure :: do
    procedure :: do_bottom
  end type type_su_mixo

contains

  subroutine initialize(self,configunit)                   
    class(type_su_mixo), intent(inout), target :: self
    integer,             intent(in)            :: configunit

    !LOCAL VARIABLES (for predation)
    integer           :: iprey
    character(len=16) :: index

    ! Register parameters in FABM
    call self%get_parameter(self%AAsyn,       'AAsyn',       'gC (gN)-1',                  ' cost for amino acid synthesis',                                              default=1.5_rk)
    call self%get_parameter(self%AEmax,       'AEmax',       'dl',                         ' maximum assimilation efficiency',                                             default=0.8_rk)
    call self%get_parameter(self%AEmin,       'AEmin',       'dl',                         ' minimum assimilation efficiency',                                             default=0.2_rk)
    call self%get_parameter(self%mBR,         'mBR',         'gC (gC-1) day-1',            ' basal respiration rate',                                                   default=0.05_rk)
    call self%get_parameter(self%mDOCpc,      'mDOCpc',      '%',                          ' % of Cfix released',                                               default=0.1_rk)
    call self%get_parameter(self%mChlCabs,    'mChlCabs',    'gChl (gC)-1',                ' absolute maximum Chl:C',                                                  default=0.03_rk)
    call self%get_parameter(self%FCabs,       'FCabs',       'gC (gC)-1',                  ' maximum feeding vacuole size',                                                   default=0.2_rk)
    call self%get_parameter(self%FCmin,       'FCmin',       'gC (gC)-1',                  ' minimum feeding vacuole size',                                                   default=0._rk)
    call self%get_parameter(self%Has,         'Has',         'dl',                         ' Hill number for digestion rate',                                                   default=2._rk)
    call self%get_parameter(self%Heq,         'Heq',         'dl',                         ' Hill number for quantity-linked AE',                                                     default=4._rk)
    call self%get_parameter(self%Hhet,        'Hhet',        'dl',                         ' Hill number for control of FCmax',                                                  default=20._rk)
    call self%get_parameter(self%Hing,        'Hing',        'dl',                         ' Hill number for ingestion control',                                                default=4._rk)
    call self%get_parameter(self%Hpbal,       'Hpbal',       'dl',                         ' Hill number for digestion link to crucial C-fixation',                                               default=4._rk)
    call self%get_parameter(self%Hpd,         'Hpd',         'dl',                         ' Hill number for digestion suppression',                                            default=10._rk)
    call self%get_parameter(self%Hrc,         'Hrc',         'dl',                         ' Hill number for regeneration',                                           default=4._rk)
    call self%get_parameter(self%mQh,         'mQh',         'dl',                         ' Hill number for control of nutrient uptake',                                                 default=4._rk)
    call self%get_parameter(self%mAKt,        'mAKt',        'ugN L-1',                    ' half-saturation for ammonium uptake',                                                 default=14._rk)
    call self%get_parameter(self%Kas,         'Kas',         'dl',                         ' half-saturation for digestion rate',                                                   default=0.5_rk)
    call self%get_parameter(self%Keq,         'Keq',         'dl',                         ' response control to ingestion quantity; 10-6 to turn off',                                                    default=0.1_rk)
    call self%get_parameter(self%Khet,        'Khet',        'dl',                         ' half-saturation for FCmax',                                                  default=1._rk)
    call self%get_parameter(self%King,        'King',        'dl',                         ' half-saturation for ingestion control',                                                default=0.2_rk)
    call self%get_parameter(self%mNKt,        'mNKt',        'ugN L-1',                    ' half-saturation for nitrate uptake',                                                 default=14._rk)
    call self%get_parameter(self%mPKu,        'mPKu',        'ugP L-1',                    ' half-saturation for phosphate uptake',                                                 default=31._rk)
    call self%get_parameter(self%KPbal,       'KPbal',       'dl',                         ' half-saturation for digestion link to crucial C-fixation',                                               default=0.1_rk)
    call self%get_parameter(self%Kpd,         'Kpd',         'dl',                         ' half-saturation for digestion supression',                                             default=1._rk)
    call self%get_parameter(self%Krc,         'Krc',         'dl',                         ' half-saturation for control of X regeneration vs re-assimilation; regeneration becomes likely as X:C tends to XCabs', default=1._rk)
    call self%get_parameter(self%Kec,         'Kec',         'dl',                         ' response control to prey quality: 10-6 to turn off',                                                    default=10._rk)
    call self%get_parameter(self%mKQN,        'mKQN',        'dl',                         ' half-saturation for N quota curve',                                                  default=10._rk)
    call self%get_parameter(self%mKxi,        'mKxi',        'dl',                         ' half-saturation for control of nutrient uptake',                                                 default=0.1_rk)
    call self%get_parameter(self%KQP,         'KQP',         'dl',                         ' half-saturation for P quota curve',                                                  default=0.1_rk)
    call self%get_parameter(self%mM,          'mM',          'dl',                         ' multiplier for Chl synthesis',                                              default=3._rk)
    call self%get_parameter(self%mMR,         'mMR',         'gC (gC)-1',                  ' metabolic respiration',                                            default=0.2_rk)
    call self%get_parameter(self%NCabs,       'NCabs',       'gN (gC)-1',                  ' absolute maximum N:C',                                                    default=0.25_rk)
    call self%get_parameter(self%NCm,         'NCm',         'gN (gC)-1',                  ' maximum N:C which could be attained in the organic form',                                                   default=0.3_rk)
    call self%get_parameter(self%mNCm,        'mNCm',        'gN (gC)-1',                  ' maximum N:C affecting growth rate',                                                   default=0.2_rk)
    call self%get_parameter(self%mNCo,        'mNCo',        'gN (gC)-1',                  ' minimum N:C',                                                    default=0.05_rk)
    call self%get_parameter(self%mPCm,        'mPCm',        'gP (gC)-1',                  ' maximum P:C affecting growth rate',                                                   default=0.02_rk)
    call self%get_parameter(self%Pbalcrit,    'Pbalcrit',    'dl',                         ' minimum critical proportion of growth supported by photosynthesis',                                         default=0._rk)
    call self%get_parameter(self%PCabs,       'PCabs',       'gP (gC)-1',                  ' absolute maximum P:C',                                                    default=0.04_rk)
    call self%get_parameter(self%PCm,         'PCm',         'gP (gC)-1',                  ' maximum P:C which could be obatined in the organic form',                                                   default=0.03_rk)
    call self%get_parameter(self%PCo,         'PCo',         'gP (gC)-1',                  ' minimum P:C',                                                    default=0.005_rk)
    call self%get_parameter(self%mAPref,      'mAPref',      'dl',                         ' relative preference for ammonium; controls surge',                                                  default=2._rk)
    call self%get_parameter(self%mNPref,      'mNPref',      'dl',                         ' relative preference for nitrate; controls surge',                                                  default=1._rk)
    call self%get_parameter(self%redco,       'redco',       'gC (gN)-1',                  ' cost of nitrate reduction to ammonium',                                               default=1.71_rk)
    call self%get_parameter(self%Smix,        'Smix',        'dl',                         'switch to mix C input; 0 if substitutional; 1 if additive mixotrophic interaction',                                            default=0)
    call self%get_parameter(self%Spd,         'Spd',         'dl',                         'switch to relate C demand to digestion; 0 if C-fixation does not affect digestion; 1 otherwise',                                 default=0)
    call self%get_parameter(self%Stype,       'Stype',       'dl',                         'switch for mixotroph type; 1 for growth; 2 for N; 3 for P; 4 for min N,P',                                                    default=1)
    call self%get_parameter(self%Svol,        'Svol',        'dl',                         'switch to share cell volume; 0 if feeding vacuole does not compromise ChlCmax',                                                default=0)
    call self%get_parameter(self%malpha,      'malpha',      '(m2/gChla)*(gC/umolphoton)', ' initial slope of PE curve',                                                  default=5e-6_rk)
    call self%get_parameter(self%aalpha,      'aalpha',      '(m2/gChla)*(gC/umolphoton)', 'algae initial slope of PE curve',                                                  default=5e-6_rk)
    call self%get_parameter(self%mbeta,       'mbeta',       'dl',                         ' control constant for nutrient uptake',                                                 default=0.05_rk)
    call self%get_parameter(self%mUm,         'mUm',         'day-1',                      ' maximum possible growth rate',                                                   default=0.7_rk)
    call self%get_parameter(self%mUphot,      'mUphot',      'day-1',                      ' maximum phototrophic growth rate',                                                   default=0.35_rk)
    call self%get_parameter(self%DChl,        'DChl',        'day-1',                      'decay rate of kleptochloroplasts',                                        default=0.693_rk*1.25_rk)
    call self%get_parameter(self%MetMult,     'MetMult',     'dl',                         ' metabolic cost multiplier to achieve required GGE for dinos',                                                  default=3._rk)
    call self%get_parameter(self%mPpref,      'mPpref',      'dl',                         'preference for DIP',                                                      default=5._rk)
    call self%get_parameter(self%iopABS,      'iopABS',      'm2/mgChl',                   'specific shortwave absorption',                                             default=0.008_rk)
    call self%get_parameter(self%iopBBS,      'iopBBS',      'm2/mgChl',                   'specific shortwave backscatter',                                            default=0._rk)
    call self%get_parameter(self%mESD,        'mESD',        'um',                         ' size in ESD',                                            default=20._rk)
    call self%get_parameter(self%pi,          'pi',          'dl',                         'number pi',                                                  default=3.14159265_rk)
    call self%get_parameter(self%wSCEB,       'wSCEB',       'm s-1',                      'root-mean-squared turbulence',                                            default=0._rk)
    call self%get_parameter(self%aSCEB,       'aSCEB',       'dl',                         'constant a for the allometric relationship of ESD and C-content', default=0.216_rk)
    call self%get_parameter(self%bSCEB,       'bSCEB',       'dl',                         'constant b for the allometric relationship of ESD and C-content', default=0.939_rk)
    call self%get_parameter(self%q10,         'q10',         '-',                          'Q_10 temperature coefficient',                                            default = 2._rk)
    call self%get_parameter(self%sd,          'sd',          'day-1',                      'specific mortality',                                                default = 0.02_rk)
    call self%get_parameter(self%Smin,        'Smin',        'um',                         'minimum prey size',                                                  default = 5._rk)
    call self%get_parameter(self%Smax,        'Smax',        'um',                         'maximum prey size',                                                  default = 20._rk)
    call self%get_parameter(self%Sopt,        'Sopt',        'um',                         'optimum prey size',                                                  default = 10._rk)
    call self%get_parameter(self%kpmax,       'kpmax',       'gC/gC/day',                  'kleptochloroplast photosynthetic rate',                             default = 0._rk)
    call self%get_parameter(self%diss_frac,   'diss_frac',   'dl',                         'dissolved fraction of detritus',                                    default = 1._rk)
    call self%get_parameter(self%sink,        'sink',        'm/day',                        'phytoplankton sinking rate',                                    default = 0._rk)

    self%dt = 86400._rk   ! Change time unit of all parameters to meet FABM default (seconds)

    ! Register internal state variables in FABM
    call self%register_state_variable(self%id_mF,      'mF',    'ugC L-1',     'food vacuole C-biomass')
    call self%register_state_variable(self%id_mFN,     'mFN',   'ugN L-1',     'food vacuole N-biomass')
    call self%register_state_variable(self%id_mFP,     'mFP',   'ugP L-1',     'food vacuole P-biomass')
    call self%register_state_variable(self%id_mFChl,   'mFChl', 'ugChl L-1',   'food vacuole Chl-biomass')
    call self%register_state_variable(self%id_mC,      'mC',    'ugC L-1',     'core  C-biomass')
    call self%register_state_variable(self%id_mN,      'mN',    'ugN L-1',     'core  N-biomass')
    call self%register_state_variable(self%id_mP,      'mP',    'ugP L-1',     'core  P-biomass')
    call self%register_state_variable(self%id_mChl,    'mChl',  'ugChl L-1',   'core  Chl-biomass') 

    ! Register variables related to the description of DIATOMS
    call self%get_parameter(self%SSi, 'SSi', '', 'switch to activate silicate limitation', default = .false.)

    if (self%SSi) then
       call self%get_parameter(self%KSi,                      'KSi',      'ugSi/L',         'half-saturation constant for silicon assimilation', default=15._rk)
       call self%get_parameter(self%mSCm,                     'mSCm',     'gSi (gC)-1',     'maximum internal Si:C',                             default=0.2_rk)
       call self%get_parameter(self%SCo,                      'SCo',      'gSi (gC)-1',     'minimum internal Si:C',                             default=0.05_rk)
       call self%get_parameter(self%SCabs,                    'SCabs',    'gSi (gC)-1',     'max absolute internal Si:C',                        default=2._rk)
       call self%get_parameter(self%betaSi,                   'betaSi',   'dl',             'constant for Si uptake control',                    default=0.4_rk)
       call self%register_state_variable(self%id_mSi,         'mSi',      'ugSi L-1',       'core  Si-biomass')
       call self%register_state_dependency(self%id_N5s,       'N5ss',     'mmolSi/m^3',     'external silicate')
       call self%register_state_dependency(self%id_bsi,       'bsi',      'mmolSi/m^2',     'benthic silicate')  
       call self%request_coupling_to_model(self%id_N5s,       'N5s',      standard_variables%total_silicate)
       call self%register_diagnostic_variable(self%id_mSC,    'mSC',      'gSi/gC',         'SiC internal quota')
       call self%register_diagnostic_variable(self%id_SCu,    'SCu',      'dl' ,            'Si limitation (0-1)')
       call self%register_diagnostic_variable(self%id_Sup,    'Sup',      'gSi/gC/d',       'Si uptake rate')
       call self%add_to_aggregate_variable(standard_variables%total_silicate,   self%id_mSi,  scale_factor=1._rk/28._rk)
    end if 

    ! Determine n prey types and its components (register also parameters related to kleptoplasty).
    call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)
    allocate(self%id_prey(self%nprey))
    allocate(self%id_preyc(self%nprey))
    allocate(self%id_preyn(self%nprey))
    allocate(self%id_preyp(self%nprey))
    allocate(self%id_preyf(self%nprey))
    allocate(self%id_preys(self%nprey))
    allocate(self%id_preyl(self%nprey))
    allocate(self%id_preychl(self%nprey))
    allocate(self%id_preyf_target(self%nprey))
    allocate(self%preyESD(self%nprey))
    allocate(self%pref(self%nprey))
    allocate(self%steal(self%nprey))
    allocate(self%id_sI(self%nprey))
    allocate(self%id_Cp(self%nprey))
    allocate(self%id_I(self%nprey))
    allocate(self%id_Ccon(self%nprey))
    allocate(self%id_Icap(self%nprey))

    do iprey=1,self%nprey
       write (index,'(i0)') iprey
       call self%get_parameter(self%preyESD(iprey),'preyESD'//trim(index),'um','size in ESD'//trim(index))
       call self%get_parameter(self%pref(iprey),'pref'//trim(index),'dl','prey preference'//trim(index))
       call self%get_parameter(self%steal(iprey),'steal'//trim(index),'','steal or not chloroplasts from prey'//trim(index))
       call self%register_dependency(self%id_preyc(iprey),'prey'//trim(index)//'c','mmol C/m^3', 'prey '//trim(index)//' carbon')
       call self%register_dependency(self%id_preyn(iprey),'prey'//trim(index)//'n','mmol N/m^3', 'prey '//trim(index)//' nitrogen')
       call self%register_dependency(self%id_preyp(iprey),'prey'//trim(index)//'p','mmol P/m^3', 'prey '//trim(index)//' phosphorus')
       call self%register_dependency(self%id_preys(iprey),'prey'//trim(index)//'s','mmol Si/m^3','prey '//trim(index)//' silicate')
       call self%register_dependency(self%id_preyl(iprey),'prey'//trim(index)//'l','mg C/m^3',   'prey '//trim(index)//' calcite')
       call self%register_dependency(self%id_preychl(iprey),'prey'//trim(index)//'chl','mg Chl/m^3',   'prey '//trim(index)//' chlorophyll')

       call self%register_diagnostic_variable(self%id_sI(iprey),   'sI'//trim(index),   'day-1',   'specific predation rate'//trim(index))
       call self%register_diagnostic_variable(self%id_Cp(iprey),   'Cp'//trim(index),   'gC/gC/d', 'capture rate'//trim(index))
       call self%register_diagnostic_variable(self%id_I(iprey),    'I'//trim(index),    'gC/gC/d', 'ingestion rate'//trim(index))
       call self%register_diagnostic_variable(self%id_Ccon(iprey), 'Ccon'//trim(index), 'dl',      'capture control'//trim(index))
       call self%register_diagnostic_variable(self%id_Icap(iprey), 'Icap'//trim(index), 'gC/gC/d', 'capture rate'//trim(index))

       call self%register_model_dependency(self%id_prey(iprey),'prey'//trim(index))
       call self%request_coupling_to_model(self%id_preyc(iprey),self%id_prey(iprey),standard_variables%total_carbon)
       call self%request_coupling_to_model(self%id_preyn(iprey),self%id_prey(iprey),standard_variables%total_nitrogen)
       call self%request_coupling_to_model(self%id_preyp(iprey),self%id_prey(iprey),standard_variables%total_phosphorus)
       call self%request_coupling_to_model(self%id_preys(iprey),self%id_prey(iprey),standard_variables%total_silicate)
       call self%request_coupling_to_model(self%id_preyl(iprey),self%id_prey(iprey),total_calcite_in_biota)
       call self%request_coupling_to_model(self%id_preychl(iprey),self%id_prey(iprey),total_chlorophyll)

       if (use_iron) then
          call self%register_dependency(self%id_preyf(iprey),'prey'//trim(index)//'f','umol Fe/m^3','prey '//trim(index)//' iron')
          call self%register_state_dependency(self%id_preyf_target(iprey),'prey'//trim(index)//'f_sink','umol Fe/m^3','sink for Fe of prey '//trim(index),required=.false.)
          call self%request_coupling_to_model(self%id_preyf(iprey),self%id_prey(iprey),standard_variables%total_iron)
       end if
    end do
 
    ! Register links to external variables (nutrient pools)
    call self%register_state_dependency(self%id_NH4,   'NH4n',   'mmolN/m^3',     'ammonium')
    call self%register_state_dependency(self%id_NO3,   'NO3n',   'mmolN/m^3',     'nitrate')
    call self%register_state_dependency(self%id_DIP,   'DIPp',   'mmolP/m^3',     'dissolved inorganic phosphorus')
    call self%register_state_dependency(self%id_DIC,   'DICc',   'mmolC/m^3',     'dissolved inorganic carbon')
    call self%register_state_dependency(self%id_sDOMC, 'sDOMCc', 'mmolC/m^3',     'semi-labile dissolved organic material-carbon')
    call self%register_state_dependency(self%id_sDOMN, 'sDOMNn', 'mmolN/m^3',     'semi-labile dissolved organic material-nitrogen')
    call self%register_state_dependency(self%id_sDOMP, 'sDOMPp', 'mmolP/m^3',     'semi-labile dissolved organic material-phosphorus')
    call self%register_state_dependency(self%id_DOC,   'DOCc',   'mmolC/m^3',     'dissolved organic carbon')
    call self%register_state_dependency(self%id_RPs,   'RPss',   'mmolSi/m^3',    'particulate organic silicate')
    call self%register_state_dependency(self%id_POC,   'POCc',   'mmolC/m^3',     'particulate organic carbon')
    call self%register_state_dependency(self%id_PON,   'PONn',   'mmolN/m^3',     'particulate organic nitrogen')
    call self%register_state_dependency(self%id_POP,   'POPp',   'mmolP/m^3',     'particulate organic phosphorus')   

    call self%register_state_dependency(self%id_bc,    'bc',     'mmolC/m^2',     'benthic carbon')
    call self%register_state_dependency(self%id_bn,    'bn',     'mmolN/m^2',     'benthic nitrogen')
    call self%register_state_dependency(self%id_bp,    'bp',     'mmolP/m^2',     'benthic phosphorus')  
 
    call self%request_coupling_to_model(self%id_NH4,   'NH4',   standard_variables%total_nitrogen)
    call self%request_coupling_to_model(self%id_NO3,   'NO3',   standard_variables%total_nitrogen)
    call self%request_coupling_to_model(self%id_DIP,   'DIP',   standard_variables%total_phosphorus)
    call self%request_coupling_to_model(self%id_DIC,   'DIC',   standard_variables%total_carbon)
    call self%request_coupling_to_model(self%id_sDOMC, 'sDOMC', standard_variables%total_carbon)
    call self%request_coupling_to_model(self%id_sDOMN, 'sDOMN', standard_variables%total_nitrogen)
    call self%request_coupling_to_model(self%id_sDOMP, 'sDOMP', standard_variables%total_phosphorus)
    call self%request_coupling_to_model(self%id_DOC,   'DOC',   standard_variables%total_carbon)
    call self%request_coupling_to_model(self%id_RPs,   'RPs',   standard_variables%total_silicate)

    if (self%diss_frac /= 1.0_rk) then
       ! Part of waste is particulate; couple to POM pools
       call self%request_coupling_to_model(self%id_POC, 'POC', standard_variables%total_carbon)
       call self%request_coupling_to_model(self%id_PON, 'PON', standard_variables%total_nitrogen)
       call self%request_coupling_to_model(self%id_POP, 'POP', standard_variables%total_phosphorus)
    else      
       ! All waste is dissolved; no need for explicitly coupled POM pools (substitute DOM)
       call self%request_coupling(self%id_POC, 'sDOMCc')
       call self%request_coupling(self%id_PON, 'sDOMNn')
       call self%request_coupling(self%id_POP, 'sDOMPp')
    end if

    if (self%sink /= 0._rk) then
       ! Sinking is active: we need a benthic POM pool to sediment into
       call self%request_coupling_to_model(self%id_bc, 'ben', standard_variables%total_carbon)
       call self%request_coupling_to_model(self%id_bn, 'ben', standard_variables%total_nitrogen)
       call self%request_coupling_to_model(self%id_bp, 'ben', standard_variables%total_phosphorus)
       if (self%SSi) call self%request_coupling_to_model(self%id_bsi, 'ben', standard_variables%total_silicate)
    else
       ! No sinking, so no need for coupled benthic pool (substitute dummies)
       call self%request_coupling(self%id_bc, 'zero_hz')
       call self%request_coupling(self%id_bn, 'zero_hz')
       call self%request_coupling(self%id_bp, 'zero_hz')
       if (self%SSi) call self%request_coupling(self%id_bsi, 'zero_hz')
    end if

    ! Register diagnostic variables
    call self%register_diagnostic_variable(self%id_mNC,    'mNC',     'gN (gC)-1',         'core NC quota')
    call self%register_diagnostic_variable(self%id_mPC,    'mPC',     'gP (gC)-1',         'core PC quota')
    call self%register_diagnostic_variable(self%id_mChlC,  'mChlC',   'gChl (gC)-1',       'core  ChlC quota')
    call self%register_diagnostic_variable(self%id_mFC,    'mFC',     'gC (gC)-1',         'feeding vacuole CC quota')
    call self%register_diagnostic_variable(self%id_mFChlC, 'mFChlC',  'gChl (gC)-1',       'feeding vacuole ChlC quota')
    call self%register_diagnostic_variable(self%id_mFNC,   'mFNC',    'gN (gC)-1',         'feeding vacuole NC quota')
    call self%register_diagnostic_variable(self%id_mFPC,   'mFPC',    'gP (gC)-1',         'feeding vacuole PC quota')
    call self%register_diagnostic_variable(self%id_UmpN,   'UmpN',    'gN/gC/day',         'maximum possible nitrogen uptake')
    call self%register_diagnostic_variable(self%id_mAVP,   'mAVP',    'gN/gC/day',         'maximum ammonium uptake')
    call self%register_diagnostic_variable(self%id_mNVP,   'mNVP',    'gN/gC/day',         'maximum nitrate uptake')
    call self%register_diagnostic_variable(self%id_mAV,    'mAV',     'gN/gC/day',         'actual ammonium uptake')
    call self%register_diagnostic_variable(self%id_mNV,    'mNV',     'gN/gC/day',         'actual nitrate uptake')
    call self%register_diagnostic_variable(self%id_mfrat,  'mfrat',   'dl',                'f-ratio')
    call self%register_diagnostic_variable(self%id_mNCu,   'mNCu',    'dl',                'normalised NC quota description')
    call self%register_diagnostic_variable(self%id_mPCu,   'mPCu',    'dl',                'normalised PC quota description')
    call self%register_diagnostic_variable(self%id_mNPCu,  'mNPCu',   'dl',                'normalised NPC quota description')
    call self%register_diagnostic_variable(self%id_mfrNC,  'mfrNC',   'dl',                'normalised feedback from NC quota')
    call self%register_diagnostic_variable(self%id_mupN,   'mupN',    'gC/gC/day',         'nitrogen uptake')
    call self%register_diagnostic_variable(self%id_mupN_mC,'mupN_mC', 'ugC/L/day',         'nitrogen uptake')
    call self%register_diagnostic_variable(self%id_mNH4up, 'mNH4up',  'ugC/L/day',         'ammonium uptake')
    call self%register_diagnostic_variable(self%id_mNO3up, 'mNO3up',  'ugC/L/day',         'nitrate uptake')
    call self%register_diagnostic_variable(self%id_mfrPC,  'mfrPC',   'dl',                'normalised feedback from PC quota')
    call self%register_diagnostic_variable(self%id_mIPu,   'mIPu',    'gP/gC/day',         'DIP uptake')
    call self%register_diagnostic_variable(self%id_mIPu_mC,'mIPu_mC', 'ugP/L/day',         'DIP uptake')
    call self%register_diagnostic_variable(self%id_PSrel,  'PSrel',   'dl',                'relative rate of photosynthesis')
    call self%register_diagnostic_variable(self%id_mPS,    'mPS',     'gC (gC)-1 (day)-1', 'constitutive-photosynthesis rate by Smith negative exponential equation')
    call self%register_diagnostic_variable(self%id_mCas,   'mCas',    'gC (gC)-1 (day)-1', 'assimilation of carbon')
    call self%register_diagnostic_variable(self%id_Rhet,   'Rhet',    'gC (gC)-1 (day)-1', 'respiration rate from heterotrophic processes')
    call self%register_diagnostic_variable(self%id_Rphot,  'Rphot',   'gC (gC)-1 (day)-1', 'phototrophic driven respiration')
    call self%register_diagnostic_variable(self%id_mCu,    'mCu',     'gC (gC)-1 (day)-1', 'net C-growth')
    call self%register_diagnostic_variable(self%id_mDgC,   'mDgC',    'gC (gC)-1 (day)-1', 'digestion rate of material held in the gut')
    call self%register_diagnostic_variable(self%id_FPS,    'FPS',     'gC (gC)-1 (day)-1', 'klepto-photosynthesis rate by Smith negative exponential equation')
    call self%register_diagnostic_variable(self%id_mPStot, 'mPStot',  'gC (gC)-1 (day)-1', 'total gross C-fix; constitutive + klepto')
    call self%register_diagnostic_variable(self%id_ImF,    'ImF',     'ugC (L)-1 (day)-1', 'total  C-ingestion')
    call self%register_diagnostic_variable(self%id_kpmax,  'kpmax',   'gC (gC)-1 (day)-1', 'max possible photosynthetic rate of kleptochloroplasts')
    call self%register_diagnostic_variable(self%id_FCmax,  'FCmax',   'gC (gC)-1',         'max possible feeding vacuoles size')
    call self%register_diagnostic_variable(self%id_Cfix,   'Cfix',    'ugC(L)-1(day)-1',   'gross photosynthesis rate')
    call self%register_diagnostic_variable(self%id_Cresp,  'Cresp',   'ugC(L)-1(day)-1',   'total respiration rate')
    call self%register_diagnostic_variable(self%id_ImFChl, 'ImFChl',  'ugChl(L)-1(day)-1', 'total chl ingestion rate')
    call self%register_diagnostic_variable(self%id_ImFk,   'ImFk',    'ugChl(L)-1(day)-1', 'chl ingestion rate from prey that provide kleptochloroplasts')
    call self%register_diagnostic_variable(self%id_photCup,'photCup', 'gC (gC)-1 (day)-1', 'C net photosynthesis')
    call self%register_diagnostic_variable(self%id_hetCup, 'hetCup',  'gC (gC)-1 (day)-1', 'C heterotrophic growth')
    call self%register_diagnostic_variable(self%id_photNup,'photNup', 'gN (gC)-1 (day)-1', 'N phototrophic growth')
    call self%register_diagnostic_variable(self%id_hetNup, 'hetNup',  'gN (gC)-1 (day)-1', 'N heterotrophic growth')
    call self%register_diagnostic_variable(self%id_photPup,'photPup', 'gP (gC)-1 (day)-1', 'P phototrophic growth')
    call self%register_diagnostic_variable(self%id_hetPup, 'hetPup',  'gP (gC)-1 (day)-1', 'P heterotrophic growth')
    call self%register_diagnostic_variable(self%id_mDgChl_mC,    'digested_Chl', 'ugChl(L)-1(day)-1', 'total chl digestion rate')
    call self%register_diagnostic_variable(self%id_Chl_decay,    'decay_Chl',    'ugChl(L)-1(day)-1', 'total kleptochloroplasts decay rate')
    call self%register_diagnostic_variable(self%id_Chl_not_used, 'Chl_not_used', 'ugChl(L)-1(day)-1', 'chl ingestion rate of non-kleptochloroplastic Chl')
    call self%register_diagnostic_variable(self%id_PAR_PB,       'PAR_PB',       'W/m2',              'light used to compute photosynthesis')
    call self%register_diagnostic_variable(self%id_mPqm,         'mPqm',         'gC/gC/day',         'potential max photosynthesis rate')
    call self%register_diagnostic_variable(self%id_Scon,         'Scon',         'dl',                'satiation control')
    call self%register_diagnostic_variable(self%id_Pbal,         'Pbal',         'dl',                'relative contribution of photosynthesis to av. growth')
    call self%register_diagnostic_variable(self%id_BalRes,       'BalRes' ,      'dl',                'relative control')
    call self%register_diagnostic_variable(self%id_PbalCon,      'PbalCon',      'dl',                'control of digestion by photosynthesis')
    call self%register_diagnostic_variable(self%id_PD,           'PD',           'dl',                'control of digestion by photosynthesis')
    call self%register_diagnostic_variable(self%id_Dmax,         'Dmax',         'gC/gC/day',         'maximum digestion rate')
    call self%register_diagnostic_variable(self%id_FCrelA,       'FCrelA',       'dl',                'relative size of feeding vacuoles')
    call self%register_diagnostic_variable(self%id_AE_op,        'AE_op',        'dl',                'actual assimilation efficiency')
    call self%register_diagnostic_variable(self%id_Urel,         'Urel',         'dl',                'relative growth rate')
    call self%register_diagnostic_variable(self%id_mRTot,        'mRTot',        'gC/gC/day',         'total respiration rate')
    call self%register_diagnostic_variable(self%id_et,           'et',           'dl',                'temperature effect')
    call self%register_diagnostic_variable(self%id_mDOCout,      'mDOCout',      'ugC/L/day',         'DOC production (a fixed fraction of photosynthesis)')
    call self%register_diagnostic_variable(self%id_mVOCL,        'mVOCL',        'ugC/L/d',           'DOC production by non-assimilated food')
    call self%register_diagnostic_variable(self%id_mVONL,        'mVONL',        'ugN/L/d',           'DON production by non-assimilated food')
    call self%register_diagnostic_variable(self%id_mVOPL,        'mVOPL',        'ugP/L/d',           'DOP production by non-assimilated food')
    call self%register_diagnostic_variable(self%id_mortC,        'mortC',        'ugC/L/d',           'DOC production due to natural mortality')
    call self%register_diagnostic_variable(self%id_mortN,        'mortN',        'ugN/L/d',           'DON production due to natural mortality')
    call self%register_diagnostic_variable(self%id_mortP,        'mortP',        'ugP/L/d',           'DOP production due to natural mortality')
    call self%register_diagnostic_variable(self%id_mregN,        'mregN',        'ugN/L/d',           'regeneration of ammonium')
    call self%register_diagnostic_variable(self%id_mregP,        'mregP',        'ugP/L/d',           'regeneration of phosphate')
    call self%register_diagnostic_variable(self%id_Preas,        'Preas',        'gP/gC/d',           'P re-assimilation')
    call self%register_diagnostic_variable(self%id_NH4reas,      'NH4reas',      'gN/gC/d',           'NH4 re-assimilation')
    call self%register_diagnostic_variable(self%id_npp,          'npp',          'ugC/L/day',         'population net primary production')
    call self%register_diagnostic_variable(self%id_FCrelV,       'FCrelV',       'dl',                'relative capacity of feeding vacuoles')
    call self%register_diagnostic_variable(self%id_stress,       'stress',       'dl',                'normalised stress controlling natural mortality')
    call self%register_diagnostic_variable(self%id_mortmC,       'mortmC',       'ugC/L/day',         'population C-mortality')

    ! Register links to external dependencies 
    call self%register_dependency(self%id_mCudep,   'mCu',    'gC (gC)-1 (day)-1', 'enables reading the previous value of the diagnostic mCu')
    call self%register_dependency(self%id_mPStotdep,'mPStot', 'gC (gC)-1 (day)-1', 'enables reading the previous value of the diagnostic mPStot')
    call self%register_dependency(self%id_avgCu,    temporal_mean(self%id_mCudep,    period=86400._rk,resolution=3600._rk,missing_value=0.0_rk))
    call self%register_dependency(self%id_avgPStot, temporal_mean(self%id_mPStotdep, period=86400._rk,resolution=3600._rk,missing_value=0.0_rk))
    call self%register_dependency(self%id_PAR, standard_variables%downwelling_photosynthetic_radiative_flux) ! local PAR
    call self%register_dependency(self%id_ETW, standard_variables%temperature)

    ! Register the contribution of variables to total N, P, C, Chl budgets
    call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_mC,  scale_factor=1._rk/12._rk)
    call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_mF,  scale_factor=1._rk/12._rk)
    call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_mN,  scale_factor=1._rk/14._rk)
    call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_mFN, scale_factor=1._rk/14._rk)
    call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_mP,  scale_factor=1._rk/31._rk)
    call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_mFP, scale_factor=1._rk/31._rk)
    call self%add_to_aggregate_variable(total_chlorophyll,                   self%id_mChl)  ! no need for scale_factor (ug/L = mg/m3)
    call self%add_to_aggregate_variable(total_chlorophyll,                   self%id_mFChl) ! no need for scale_factor (ug/L = mg/m3)

    ! Register the contribution to fluxes
    call self%add_to_aggregate_variable(photosynthesis_rate,                 self%id_Cfix)  ! no need to convert units (mgC/m3/d = ugC/L/d)

    ! Register contribution to light extinction
    call self%add_to_aggregate_variable(particulate_organic_absorption_coefficient,  self%id_mChl,  scale_factor=self%iopABS)      
    call self%add_to_aggregate_variable(particulate_organic_backscatter_coefficient, self%id_mChl,  scale_factor=self%iopBBS)   
    call self%add_to_aggregate_variable(particulate_organic_absorption_coefficient,  self%id_mFChl, scale_factor=self%iopABS)     
    call self%add_to_aggregate_variable(particulate_organic_backscatter_coefficient, self%id_mFChl, scale_factor=self%iopBBS)  

  end subroutine initialize

  subroutine do(self, _ARGUMENTS_DO_)
      class (type_su_mixo), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ETW, et, mUm, mUphot                                              ! temperature effect
      real(rk) :: mF, mFN, mFP, mFChl, mC, mN, mP, mChl, mSi                        ! internal state variables
      real(rk) :: NH4, NO3, DIP, DIC, sDOMC, sDOMN, sDOMP, DOC, RPs, N5s            ! external state variables
      real(rk) :: POC, PON, POP                                                     ! external state variables
      real(rk) :: avgCu, avgPStot                                                   ! external dependencies
      integer  :: iprey, istate                                                     ! PREDATION: n prey
      real(rk),dimension(self%nprey) :: preycP,preypP,preynP,preysP,preylP,preychlP ! PREDATION: nutrient components within prey i
      real(rk),dimension(self%nprey) :: prey_radius, preyCcell, NAprey, SMprey      ! PREDATION: SCEB function
      real(rk) :: m_radius, mCcell, SMm                                             ! PREDATION: SCEB function
      real(rk),dimension(self%nprey) :: ENCs, ENCd, Cp, Ccon, Icap, I, sI           ! PREDATION: SCEB function and ingestion rates on prey i
      real(rk) :: preyP, Chl_not_used                                               ! PREDATION
      real(rk) :: mNC, mPC, mChlC, mFChlC                                           ! protist quotas
      real(rk) :: UmpN, mAVP, mNVP, mAV, mNV, mfrat                                 ! nitrogen transport
      real(rk) :: mNCu, mPCu, mNPCu                                                 ! N and P limitation
      real(rk) :: mfrNC, mupN, mupN_mC, mNH4up, mNO3up, mfrPC, mIPu, mIPu_mC        ! nitrogen uptake
      real(rk) :: mSC, VS, SCu, mfrSi, Sup, mortS                                   ! for DIATOMS only
      real(rk) :: BRop, mPqm, PAR, mPS                                              ! constitutive photosynthesis
      real(rk) :: kpmax, FPS, mPStot, PSrel, Cfix                                   ! non-constitutive photosynthesis
      real(rk) :: Urel, FCmax, mFC, FCrelV, Scon                                    ! INGESTION
      real(rk) :: sumCp, ImF, ImFN, ImFP, ImFChl, ImFk                              ! INGESTION
      real(rk) :: mFNC, mFPC, RelNC, RelPC, MINup, AEqual, AEquan, mAE, AE_op       ! Digestion
      real(rk) :: Dmax, PD, FCrelA, Pbal, BalRes, PbalCon                           ! Digestion 
      real(rk) :: mDgC, mDgC_mC, mDgN_mC, mDgP_mC, mDgChl_mC, Chl_decay             ! Digestion
      real(rk) :: mCas, AgC_mC, IncN, IncN_mC, IncP, IncP_mC                        ! Assimilation
      real(rk) :: XSC, mBRi, BRb, Rhet                                              ! Heterotrophic respiration
      real(rk) :: mNrep, RegN, NH4reas, Prep, RegP, Preas                           ! Re-assi regenerated N/P
      real(rk) :: Rphot, mRTot, Cresp, mICout, mCu                                  ! Respiration / growth
      real(rk) :: Cin, ChlCm, dChlC, dChlC_mC                                       ! Chl synthesis
      real(rk) :: resN, resP, resN_mC, resP_mC                                      ! N and P regeneration
      real(rk) :: mfFNC, mfFPC, mVONtr, mVOPtr, mVOCt, mVONt, mVOPt, mregN, mregP   ! Non-assimilated food
      real(rk) :: mICup, mDOCout, mVOCL, mVONL, mVOPL                               ! DIC usage and voided sC,sN,sP
      real(rk) :: photCup, hetCup, photNup, hetNup, photPup, hetPup, PGGE, HGGE     ! photo vs hetero growth 
      real(rk) :: npp, total_steal
      real(rk) :: stress, mortmC, mortmN, mortmP, mortmFC, mortmFN, mortmFP, mortmChl, mortmFChl, mortC, mortN, mortP, DChl, sd

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

      ! Retrieve internal state variables values
      _GET_(self%id_mF,    mF)
      _GET_(self%id_mFN,   mFN)
      _GET_(self%id_mFP,   mFP)
      _GET_(self%id_mFChl, mFChl)
      _GET_(self%id_mC,    mC)
      _GET_(self%id_mN,    mN)
      _GET_(self%id_mP,    mP)
      _GET_(self%id_mChl,  mChl)

      ! Retrieve external state variables values
      _GET_(self%id_NH4,   NH4)
      _GET_(self%id_NO3,   NO3)
      _GET_(self%id_DIP,   DIP)
      _GET_(self%id_DIC,   DIC)
      _GET_(self%id_DOC,   DOC)
      _GET_(self%id_sDOMC, sDOMC)
      _GET_(self%id_sDOMN, sDOMN)
      _GET_(self%id_sDOMP, sDOMP)
      _GET_(self%id_RPs,   RPs)
      _GET_(self%id_POC,   POC)
      _GET_(self%id_PON,   PON)
      _GET_(self%id_POP,   POP)

      ! Nutrient pools were returned in mmolX/m^3; convert to ugX/L (units used in this model)
      NH4 = NH4*14._rk
      NO3 = NO3*14._rk
      DIP = DIP*31._rk
      DIC = DIC*12._rk
      DOC = DOC*12._rk
      sDOMC = sDOMC*12._rk
      sDOMN = sDOMN*14._rk
      sDOMP = sDOMP*31._rk
      RPs = RPs*28._rk
      POC = POC*12._rk
      PON = PON*14._rk
      POP = POP*31._rk

      ! Get prey concentrations
      do iprey=1,self%nprey
         _GET_(self%id_preyc(iprey), preycP(iprey))
         _GET_(self%id_preyn(iprey), preynP(iprey))
         _GET_(self%id_preyp(iprey), preypP(iprey))
         _GET_(self%id_preys(iprey), preysP(iprey))
         _GET_(self%id_preyl(iprey), preylP(iprey))
         _GET_(self%id_preychl(iprey), preychlP(iprey))
      end do

      ! Prey was returned in mmol (due to units of e.g., standard_variables%total_carbon); convert to ugX/L (except prey chl; units are correct)
      preycP = preycP*12._rk
      preynP = preynP*14._rk
      preypP = preypP*31._rk
      preysP = preysP*28._rk

      ! Retrieve external dependencies
      _GET_(self%id_mCudep,    mCu)
      _GET_(self%id_mPStotdep, mPStot)
      _GET_(self%id_avgCu,     avgCu)
      _GET_(self%id_avgPStot,  avgPStot)
      _GET_(self%id_PAR,       PAR)
      _GET_(self%id_ETW,       ETW)

!---------Calculations:---------

         ! Temperature effect:
         et = self%q10**((ETW-10._rk)/10._rk) - self%q10**((ETW-32._rk)/3._rk)
         
         ! Temperature affects the following rates:
         mUm = self%mUm * et
         mUphot = self%mUphot * et 
         kpmax = self%kpmax * et
         sd = self%sd * et
         DChl = self%DChl * et

!------------------------------------
!----Silicon dynamics for diatoms----
!------------------------------------

      if (self%SSi) then

         _GET_(self%id_mSi,   mSi) ! Retrieve internal state variables values
         _GET_(self%id_N5s,   N5s) ! Retrieve external state variables values
         N5s = N5s*28._rk          ! values returned in mmolSi/m3; convert it to ugSi/L

         ! Internal SiC quota (gSi/gC):
         mSC = mSi/mC

         ! Potential silicon uptake rate (gSi/gC/day):
         if (N5s > 0._rk) then
            VS = mUphot * self%mSCm * N5s/(N5s + self%KSi)
         else
            VS = 0._rk
         end if

         ! Silicon limitation factor (dl):
         if (mSC > self%SCo) then
            if (VS >= mUphot*self%SCo) then
               SCu = 1._rk 
            else
               SCu = VS/(mUphot*self%SCo) 
            end if
         else
            SCu = 0._rk
         end if

         ! Normalized feedback to compute the actual silicon uptake Sup (dl):
         if (N5s > 0._rk) then
            mfrSi = N5s/(N5s + self%KSi)*(1._rk - mSC/self%SCabs)**self%mQh/((1._rk - mSC/self%SCabs)**self%mQh + self%mKxi)
         else
            mfrSi = 0._rk
         end if

      end if  

!-------------------------------------
!-----------Nutrient uptake-----------
!-------------------------------------

      ! Nitrogen transport
      !  maximum N-source uptake potential (gN/gC/day):
      UmpN = mUphot*self%mNCm

      !  maximum N-source uptake potential (gN/gC/day):
      if (NH4 > 0._rk) then 
         mAVP = UmpN * self%mApref * NH4/(NH4 + self%mAkt)
      else 
         mAVP = 0._rk
      end if

      !  potential relative NO3 transport (gN/gC/day):
      if (NO3 > 0._rk) then 
         mNVP = UmpN * self%mNpref * NO3/(NO3 + self%mNkt)
      else 
         mNVP = 0._rk
      end if 

      !  actual NH4 uptake (gN/gC/day):
      if (UmpN < mAVP) then
         mAV = UmpN
      else
         mAV = mAVP
      end if 

      !  actual NO3 transport (gN/gC/day):
      if (mAVP < UmpN) then
         if (mAVP + mNVP < UmpN) then 
            mNV = mNVP
         else
            mNV = UmpN - mAVP
         end if
      else
         mNV = 0._rk 
      end if   

      !  f-ratio (dl):
      mfrat = mNV / (mNV + mAV + 1e-12_rk)

      ! NC internal quota (gN/gC):
      mNC = mN/mC

      ! PC internal quota (gP/gC):
      mPC = mP/mC

      !  normalised NC quota description (dl):
      if (mNC <= self%mNCm) then
         if (mNC > self%mNCo) then
            mNCu = (1._rk + self%mKQN) * (mNC - self%mNCo)/((mNC - self%mNCo) + self%mKQN * (self%mNCm - self%mNCo))
         else
            mNCu = 0._rk
         end if 
      else
      mNCu = 1._rk
      end if 

      !  normalised PC quota description (dl): 
      if (mPC <= self%mPCm) then
         if (mPC > self%PCo) then
            mPCu = (1._rk + self%KQP) * (mPC - self%PCo)/((mPC - self%PCo) + self%KQP * (self%mPCm - self%PCo))
         else
            mPCu = 0._rk
         end if 
      else
      mPCu = 1._rk
      end if

      !  threshold selection of phototrophic growth control (by N or P or Si (if diatoms) ) (dl):
      if (self%SSi) then
         mNPCu = MIN(1._rk, mNCu, mPCu, SCu)
      else 
         mNPCu = MIN(1._rk, mNCu, mPCu)
      end if

      !  normalised feedback from NC quota (dl)
      mfrNC = (1._rk + self%mKxi**self%mQh)*(1._rk - mNC/self%NCabs)**self%mQh/((1._rk - mNC/self%NCabs)**self%mQh + self%mKxi**self%mQh)

      !  total DIN uptake (gN/gC/day):
      if (mNC < self%NCabs - 0.002_rk) then 
         if (mNCu > mNPCu) then
            mupN = (mAV + mNV)*mNPCu**self%mbeta*mfrNC
         else 
            mupN = (mAV + mNV)*mfrNC
         end if
      else
         mupN = 0._rk
      end if

      !  population total DIN uptake (ugN/L/day):
      mupN_mC = mupN*mC

      !  population total NH4 uptake (ugN/L/day):
      mNH4up = mupN_mC*(1._rk - mfrat)

      !  population total NO3 uptake (ugN/L/day):
      mNO3up = mupN_mC*mfrat

      !  normalised feedback from PC quota (dl)
      mfrPC = DIP/(DIP + self%mPKu)*(1._rk + self%mKxi**self%mQh)*(1._rk - mPC/self%PCabs)**self%mQh/((1._rk - mPC/self%PCabs)**self%mQh + self%mKxi**self%mQh)

      !  total DIP uptake (gP/gC/day):
      if (DIP > 0._rk .and. mPC < self%PCabs) then 
         if (mPCu > mNPCu) then
            mIPu = mUphot*self%mPpref*self%mPCm*mNPCu**self%mbeta*mfrPC
         else 
            mIPu = mUphot*self%mPpref*self%mPCm*mfrPC
         end if
      else
         mIPu = 0._rk
      end if

      !  population total DIP uptake (ugP/L/day):
      mIPu_mC = mIPu*mC   

!-----------------------------------
!----Constitutive Photosynthesis----
!-----------------------------------

      ! Basal respiration rate, including term to halt respiration at high NC (gC/gC/d):
      if (mNC < self%NCabs) then
         BRop = mUm*self%mBR*1.01_rk*((self%NCabs - mNC)/(self%NCabs - self%mNCo))/((self%NCabs - mNC)/(self%NCabs - self%mNCo) + 0.01_rk)
      else
         BRop = 0._rk
      end if

      !  Constitutive potential maximum photosynthesis rate (gC/gC/day):
      mPqm = (mUphot + BRop + self%mNCm*mUphot*(self%redco + self%AAsyn*self%MetMult))*mNPCu + 1e-6_rk

      !  ChlC quota (gChl/gC)
      mChlC = mChl/mC

      ! Constitutive photosynthesis rate (gC/gC/day):
      ! compute according to the Smith negative exponential equation:
      ! PAR (W/m2) is multiplied by 4.57 to convert it to PFD (photon flux density: umolphoton/m2/s)
      ! PFD is then multiplied by 24*60*60 to get it in units of per day
      mPS = mPqm * (1._rk - exp(-self%malpha * mChlC * PAR * 4.57_rk *24._rk*60._rk*60._rk/mPqm))

!------------------------------------
!-------------Ingestion--------------
!------------------------------------

      !  relative growth rate (dl):
      if (self%Stype == 1) then 
         Urel = MIN(1._rk, avgCu/mUm)
      else if (self%Stype == 2) then
         Urel = mNCu
      else if (self%Stype == 3) then
         Urel = mPCu
      else
         Urel = MIN(mNCu, mPCu)
      end if

      !  maximum operational gut size in response to C, N, and P stress (gC/gC):
      if (self%Smix == 0) then
         FCmax = self%FCmin + (1._rk - (1._rk + self%Khet**self%Hhet)*Urel**self%Hhet/(Urel**self%Hhet + self%Khet**self%Hhet))*(self%FCabs - self%FCmin)
      else 
         FCmax = self%FCabs
      end if

      ! Carbon in the feeding vacuole : carbon in the  core (gC/gC):
      mFC = mF/mC

      !  relative capacity of gut (used to halt ingestion as the gut fills) (dl):
      if (mFC > 0._rk) then
         FCrelV = MIN(1._rk, mFC/FCmax)
      else 
         FCrelV = 0._rk
      end if

      ! Normalized satiation control (dl):
      Scon = (1._rk + self%King**self%Hing)*(1._rk - FCrelV)**self%Hing/((1._rk - FCrelV)**self%Hing + self%King**self%Hing)

!---- Implementing allometry and prey selection (Flynn and Mitra 2016 Frontiers in Marine Science; Flynn, 2018 Book)----!

      ! Predator radius and radius of prey i (m):
      m_radius =    self%mESD/2._rk/10e6_rk
      prey_radius = self%preyESD/2._rk/10e6_rk

      !  carbon content per cell and carbon content per cell of prey i (ng C/cell):
      mCcell =    self%aSCEB/1000._rk * ((4._rk/3._rk * self%pi * self%mESD/2._rk)**3._rk)**self%bSCEB
      preyCcell = self%aSCEB/1000._rk * ((4._rk/3._rk * self%pi * self%preyESD/2._rk)**3._rk)**self%bSCEB
  
      ! Numeric prey i abundance (cells/m3):
      NAprey = preycP*10e6_rk /preyCcell

      ! Predator (SMm) and prey i (SMprey) speed of motion (m/s):
      SMm =    10e-6_rk * (38.542_rk * et * self%mESD**0.5424_rk) 

      do iprey=1,self%nprey
         if (self%preyESD(iprey) - int(self%preyESD(iprey)) == 0._rk) then
            SMprey(iprey) = 10e-6_rk * (38.542_rk * et * self%preyESD(iprey)**0.5424_rk)
         else
            SMprey(iprey) = 0._rk
         end if
      end do

      ! Encounter rate by the predator and prey i (prey/predator/s):
      ENCs = self%pi * (prey_radius + m_radius)**2._rk * NAprey * (SMprey**2._rk + 3._rk * SMm**2._rk + 4._rk * self%wSCEB**2._rk) * (SMm**2._rk + self%wSCEB**2._rk)**(-0.5_rk) * 3**(-1._rk)
      
      ! Encounter rate per day (prey/predator/day):
      ENCd = ENCs * 60._rk * 60._rk * 24._rk  

      ! Potential ingestion rate of prey i assuming that all encountered prey are actually ingested (gCprey/gCmixo/day):
      Cp = ENCd * (preyCcell/mCcell)

      ! Normalised capture control of prey i according to prey size preference (dl):
      do iprey=1,self%nprey      
         if (self%preyESD(iprey) > self%Smin .and. self%preyESD(iprey) < self%Smax) then
            if (self%preyESD(iprey) < self%Sopt) then
               Ccon(iprey) = (self%preyESD(iprey) - self%Smin)/(self%Sopt - self%Smin)
            else
               Ccon(iprey) = (self%Smax - self%preyESD(iprey))/(self%Smax - self%Sopt)
            end if
         else
            Ccon(iprey) = 0._rk
         end if
      end do

      ! Actual ingestion rate of prey i considering capture control and prey type preference (gCprey/gCmixo/day) 
      Icap = Cp * Ccon * self%pref

      ! Actual ingestion rate of prey i considering capture and satiation controls and prey type preference (I; gCprey/gCmixo/day) 
      I = Cp * Ccon * Scon * self%pref

      ! Specific ingestion rate on prey i (sI; day-1):
      do iprey=1,self%nprey
         if (preycP(iprey) > 0._rk) then
               sI(iprey) = I(iprey) * mC / preycP(iprey)
         else 
            sI(iprey) = 0._rk
         end if
      end do

      ! Total ingestion rate in ug X prey/L/day:
      ImF  =   sum(sI * preycP)
      ImFN =   sum(sI * preynP)
      ImFP =   sum(sI * preypP)
      ImFChl = sum(sI * preychlP)              ! total amount of Chl ingested
      ImFk   = sum(sI * preychlP * self%steal) ! chl ingested from prey that provide kleptochloroplasts
      Chl_not_used = ImFChl - ImFk             ! chl from prey that do not provide kleptochloroplasts 

      ! Apply specific predation rates to all state variables of every prey.
      do iprey=1,self%nprey
         do istate=1,size(self%id_prey(iprey)%state)
            _GET_(self%id_prey(iprey)%state(istate),preyP)
            _ADD_SOURCE_(self%id_prey(iprey)%state(istate),-sI(iprey)*preyP)
         end do
      end do 

!---------------------------------------
!----Non-constitutive Photosynthesis----
!---------------------------------------

      ! Chl in the feeding vacuole : carbon in the core (gChl/gC):
      if (mC > 0._rk) then
         mFChlC = mFChl/mC
      else
         mFChlC = 0._rk
      end if

      ! Photosynthesis rate associated with kleptochloroplasts (gC/gC/day):
      ! compute according to the Smith negative exponential equation:
      ! PAR (W/m2) is multiplied by 4.57 to convert it to PFD (photon flux density: umolphoton/m2/s)
      ! PFD is then multiplied by 24*60*60 to get it in units of per day
      if (kpmax == 0._rk .or. mNPCu == 0._rk) then
         FPS = 0._rk
      else
         FPS = kpmax*mNPCu * (1._rk - exp(-self%aalpha * PAR*4.57_rk *24._rk*60._rk*60._rk * mFChlC/kpmax * mNPCu))
      end if
 
      !  total gross C-fix, including constitutive + klepto (gC/gC/day):
      mPStot = mPS + FPS

      ! C-fixation rate by  population (ugC/L/day):
      Cfix = mPStot*mC + 1e-60_rk

      !  relative rate of photosynthesis (dl):
      PSrel = MIN(1._rk, mPStot/mPqm)

!--------------------------------------
!------Digestion and assimilation------
!--------------------------------------

      ! Nitrogen in the feeding vacuole : carbon in the  core (gC/gC):
      mFNC = mFN/mC

      ! Phosphorous in the feeding vacuole : carbon in the  core (gC/gC):
      mFPC = mFP/mC

      !  relative prey N to predator N (dl)
      RelNC = (mFNC/mFC)/self%mNCm

      !  relative prey P to predator P (dl)
      RelPC = (mFPC/mFC)/self%mPCm

      ! Minimal threshold control to select the release of N/P for differences in SXC and ZXC, ratio (dl)
      MINup = MIN(RelNC, RelPC, 1._rk)

      ! Efficiency parameter for assimilation/digestion, RH function linked to prey QUALITY giving a value between 0 & 1 (i.e. quotient). AE can never reach 1; (dl):
      AEqual = self%AEmin + (self%AEmax - self%AEmin)*(1._rk + self%Kec)*MINup/(MINup + self%Kec)

      ! AE changing in response to changes in relative growth rate; decreases with increases growth rate; i.e. linked to QUANTITY. Turn this off by setting KecQ=1e-6 (dl):
      if (self%Keq > 1e-6_rk) then
         AEquan = self%AEmin + (self%AEmax - self%AEmin)*(1._rk + self%Keq**self%Heq)*(1._rk - Urel)**self%Heq/((1._rk - Urel)**self%Heq + self%Keq**self%Heq)
      else
         AEquan = self%AEmax
      end if

      ! Actual AE as the minimum of quality and quantity linked AE (dl):
      mAE = MIN(AEqual, AEquan)

      ! Operational AE (dl):
      AE_op = mAE*MINup + 1e-6_rk

      !  digestion rate required to support Um with the current food quality AE_FQ (gC/gC/d):
      Dmax = (mUm + BRop)/(AE_op*(1._rk - self%mMR))

      !  control on digestion by PSrel; if Spd switch is enabled (Spd = 1) then digestion of captured prey is linked to PSrel - if sufficient C is flowing from PS then decrease digestion (dl): 
      PD = 1._rk - self%Spd*((1._rk + self%Kpd**self%Hpd)*PSrel**self%Hpd/(PSrel**self%Hpd + self%Kpd**self%Hpd))

      !  size of FC relative to the maximum possible FCabs (availability of material in the gut for digestion) (dl):
      if (mFC > 0._rk) then
         FCrelA = mFC/self%FCabs
      else
         FCrelA = 0._rk
      end if

      !  balance of total photosynthesis to growth rate (dl):
      if (avgCu > 0._rk) then
         Pbal = MIN(1._rk, avgPStot/(avgCu + 1e-24_rk))
      else 
         Pbal = 0._rk
      end if

      !  hetero vs photo balance (dl):
      if (Pbal > self%Pbalcrit) then
         BalRes = (Pbal - self%Pbalcrit)/(1._rk - self%Pbalcrit)
      else 
         BalRes = 0._rk
      end if

      !  control of digestion from material in the gut by the relative rate of photosynthesis (dl):
      if (self%Pbalcrit == 0._rk) then
         PbalCon = 1._rk
      else
         PbalCon = (1._rk + self%KPbal**self%Hpbal)*BalRes**self%Hpbal/(BalRes**self%Hpbal + self%KPbal**self%Hpbal)
      end if

      ! Digestion of material held in the gut, as a sigmoidal function of gut contents relative to the value of FCrelA (gC/gC/d):  
      if (mFC > 10e-201_rk) then
         mDgC = PbalCon*Dmax*(1._rk + self%Kas**self%Has)*FCrelA**self%Has/(FCrelA**self%Has + self%Kas**self%Has)*PD
      else
         mDgC = 0._rk
      end if

      ! Digestion by  population (ugX prey/L/day):
      mDgC_mC   = mDgC*mC
      mDgN_mC   = (mDgC/mFC)*mFN
      mDgP_mC   = (mDgC/mFC)*mFP
      mDgChl_mC = (mDgC/mFC)*mFChl

      ! Dilution of kleptochloroplasts
      Chl_decay = DChl*mFChl

      ! Assimilation of carbon (gX/gC/d and ugX/L/day):
      mCas   =  AE_op*mDgC
      AgC_mC =  mCas*mC

      if (mCas > 0._rk) then
         IncN =    mCas*self%mNCm
      else
         IncN =    0._rk
      end if

      IncN_mC = IncN*mC

      if (mCas > 0._rk) then
         IncP =    mCas*self%mPCm
      else
         IncP =    0._rk
      end if

      IncP_mC = IncP*mC

!---------------------------------------
!-----Respiration and Regeneration------
!---------------------------------------

      !  supply of C for support of respiration but not used for protoplasmic purposes (C/C/d):
      XSC = AEqual*mDgC*(1._rk - MINup)

      ! Support of BR using XSC (gC/gC/d):
      if (BRop <= XSC) then
         mBRi = BRop
      else
         mBRi = XSC
      end if

      ! Basal resp from predator body, this occurs when XSC cannot support BR (gC/gC/d):
      BRb = BRop - mBRi

      ! Respiration rate from heterotrophic processes (gC/gC/d):
      Rhet = BRb + mCas*self%mMR

      !  repressor of regenerated N reassimilation; decreases as NC approaches NCabs (dl):
      if (mNC > self%mNCm) then
         mNrep = 1._rk - (self%NCabs - mNC)/(self%NCabs - self%mNCm)
      else
         mNrep = 0._rk
      end if  

      !  control of N regeneration (dl):
      RegN = (1._rk + self%Krc**self%Hrc)*mNrep**self%Hrc/(mNrep**self%Hrc + self%Krc**self%Hrc)

      !  NH4 re-assimilated (gN/gC/d):
      NH4reas = Rhet*self%mNCm*(1._rk - RegN)

      !  repressor of regenerated P reassimilation; decreases as PC approaches PCabs (dl):   
      if (mPC > self%mPCm) then
         Prep = 1._rk - (self%PCabs - mPC)/(self%PCabs - self%mPCm)
      else
         Prep = 0._rk
      end if  

      !  control of P regeneration (dl):
      RegP = (1._rk + self%Krc**self%Hrc)*Prep**self%Hrc/(Prep**self%Hrc + self%Krc**self%Hrc)   

      !  P re-assimilated (gP/gC/d):
      Preas =  Rhet*self%mPCm*(1._rk - RegP)

      !  phototrophic driven respiration (gC/gC/d):
      Rphot = self%redco*mupN*mfrat+(mupN+NH4reas)*self%AAsyn*self%MetMult

      !  total respiration (gC/gC/d):
      mRTot = Rhet + Rphot

      ! Total respiration loss from  population (ugC/L/d):
      Cresp = mRTot*mC

      ! Total DIC out by  population (ugC/L/d):
      mICout = mC*(mRTot + mBRi)

      !  net C-growth (gC/gC/d):
      mCu = mCas + mPStot - mRTot

!---------------------------------
!----Non-predatory mortality------
!---------------------------------

      ! Normalised stress (0-1; dl):
      if (avgCu > 0._rk) then
         if (avgCu > mUm) then
             stress = 0._rk
         else
             stress = 1._rk - avgCu/mUm
         end if
      else
      stress = 0._rk
      end if

      ! Linear mortality term
      mortmC  = stress * sd *  mC   ! (ugC/L/day)
      mortmN  = stress * sd *  mN   ! (ugN/L/day)
      mortmP  = stress * sd *  mP   ! (ugP/L/day)
      mortmFC = stress * sd * mF    ! (ugC/L/day)
      mortmFN = stress * sd * mFN   ! (ugN/L/day)
      mortmFP = stress * sd * mFP   ! (ugP/L/day)
      mortmChl = stress * sd * mChl   ! (ugChl/L/day)
      mortmFChl = stress * sd * mFChl ! (ugChl/L/day)

      mortC = mortmC + mortmFC
      mortN = mortmN + mortmFN
      mortP = mortmP + mortmFP

!-------------------------------------
!-----Silicon dynamics for diatoms----
!-------------------------------------

      if (self%SSi) then

         ! Actual silicon uptake (gSi/gC/day):
         if (mSC < self%SCabs - 0.01_rk) then
            if ((SCu > mNPCu) .or. (SCu == 1._rk .and. mNPCu == 1._rk)) then      ! reads as no Si limitation
               if (mCu > 0._rk) then
                  Sup = mUphot * self%mSCm * (mCu/mUphot)**self%betaSi * mfrSi
               else          
                  Sup = 0._rk
               end if
            else                                                                  ! reads as under Si limitation
               Sup = mUphot * self%mSCm * mfrSi
            end if
         else
            Sup = 0._rk
         end if

         ! Mortality associated to Si
         mortS  = stress * sd *  mSi   ! (mgSi/m3/day)
         
         ! Set ODEs: 
         _ADD_SOURCE_(self%id_mSi, + Sup*mC - mortS)
         _ADD_SOURCE_(self%id_N5s, - (Sup*mC)/28._rk)

         ! Set Diagnostics:
         _SET_DIAGNOSTIC_(self%id_mSC, mSC)
         _SET_DIAGNOSTIC_(self%id_SCu, SCu)
         _SET_DIAGNOSTIC_(self%id_Sup, Sup)

      end if

!---------------------------------------
!---------Chlorophyll synthesis---------
!---------------------------------------

      !  total C input to control the demand for Chl; if Smix=0 then this is substitutional, with C from ingestion compensating for a shortfall in C-fixation. If Smix=1 then only the contribution by Cfix is considered (C/C/d):
      if (self%Smix == 0) then
         Cin = mPStot + mCas
      else
         Cin = mPStot
      end if

      ! Maximum Chl:C (gChl/gC); the value of 12 is C:Chl for the volume of photosystems:
      if (self%Svol == 1) then
         ChlCm = self%mChlCabs - ((mFC/(1._rk + mFC))/12._rk)
      else
         ChlCm = self%mChlCabs
      end if

      !  change in ChlC with synthesis and degradation with C growth (g Chl g-1C d-1):
      if (mChlC < ChlCm) then
         dChlC = ChlCm*mNPCu*self%mM*mUphot*(1._rk - MIN(1._rk, Cin/mPqm))**0.5_rk*(1._rk + 0.05_rk)*(1._rk - mChlC/ChlCm)/(1._rk - mChlC/ChlCm + 0.05_rk) - mChlC*((1._rk - mNCu)*mUphot)
      else
         dChlC = 0._rk
      end if

      !  change in ChlC with synthesis and degradation with C growth (ug Chl/L/day):
      dChlC_mC = dChlC*mC

!--------------------------------
!------------Voiding-------------
!--------------------------------

      !  stoichiometric-linked X respiration with heterotrophic respiration, assuming that the base X:C of core biological material accords with XCm  (gX/gC/d and ugX/L/d):
      resN = Rhet*self%mNCm*RegN
      resP = self%mPCm*Rhet*RegP
      resN_mC = resN*mC
      resP_mC = resP*mC

      ! XC of material in the feeding vacuole (gX in F : gC in mixo) / (gC in F : gC in mixo) which means gX in F / gC in F:
      mfFNC = mFNC/mFC
      mfFPC = mFPC/mFC

      !  rate of X voiding as balance of that ingested versus that required to maintain  X:C  (gX/gC/d), which means: all that was digested - what was assimilated = non-assimilated food (in units of ugX/L/d this could be obtained simply by: mDgX_mC - IncX_mC):
      mVONtr = mDgC*mfFNC - mCas*self%mNCm
      mVOPtr = mDgC*mfFPC - mCas*self%mPCm

      !  non-assimilated food of C as DOC (gC/gC/d):
      mVOCt = mDgC - mCas - (mBRi + 100e-102_rk)

      !  rate of removal of excess nitrogen/phosphorus to maintain X:C constant (gX/gC/d):
      if (mVONtr/mVOCt > self%NCm) then
         mVONt = mVOCt*self%NCm + 1e-9_rk
      else
         mVONt = mVONtr + 1e-9_rk
      end if

      if (mVOPtr/mVOCt > self%PCm) then
         mVOPt = mVOCt*self%PCm + 1e-9_rk
      else
         mVOPt = mVOPtr + 1e-9_rk
      end if

      !  regeneration rate of N as NH4 and P as DIP (gX/gC/d) * mC = (ugX/L/d):
      mregN = (resN + (mVONtr - mVONt))*mC
      mregP = (resP + (mVOPtr - mVOPt))*mC

      ! Total  DIC usage (ugC/L/day):
      mICup = mPStot*mC
      mDOCout = mICup*self%mDOCpc
      
      !  rate of removal of excess carbon to maintain the X:C constant (ugC/L/d):
      mVOCL = mVOCt*mC

      !  loss of N as PON (ugN/L/d):
      mVONL = mVONt*mC

      !  loss of P as POP (ugP/L/d):
      mVOPL = mVOPt*mC

      !  net photosynthesis rate (photCup) and heterotrophic C-assimilation (hetCup)(gC/gC/d):
      photCup = mPStot - Rphot
      hetCup  = mCas - Rhet
      !  total phototrophic (photNup) and heterotrophic (hetNup) N assimilation (gN/gC/d):
      photNup = mupN + NH4reas
      hetNup  = IncN - (resN + NH4reas)
      !  total phototrophic (photPup) and heterotrophic (hetPup) P assimilation (gP/gC/d):
      photPup = mIPu + Preas
      hetPup  = IncP - (resP + Preas)

      ! Net primary productivity (NPP) considering all C-fixation and all C-respired by the population (ugC/L/day):
      npp = Cfix - Cresp

      if (use_iron) then
         ! Iron dynamics:
         ! Following Vichi et al., 2007 it is assumed that the iron fraction of the ingested phytoplankton
         ! is egested as particulate detritus (Luca)
         do iprey=1,self%nprey
            _GET_(self%id_preyf(iprey),preyP)
            if (preyP/=0.0_rk) _ADD_SOURCE_(self%id_preyf_target(iprey),sI(iprey)*preyP)
         end do
      end if

!-----------------------------
!---Set ODEs and diagnostics--
!-----------------------------

      _ADD_SOURCE_(self%id_mC,    + Cfix     + AgC_mC   - Cresp   - mortmC)
      _ADD_SOURCE_(self%id_mN,    + mupN_mC  + IncN_mC  - resN_mC - mortmN)
      _ADD_SOURCE_(self%id_mP,    + mIPu_mC  + IncP_mC  - resP_mC - mortmP)
      _ADD_SOURCE_(self%id_mChl,  + dChlC_mC - mortmChl)
      _ADD_SOURCE_(self%id_mF,    + ImF      - mDgC_mC - mortmFC)
      _ADD_SOURCE_(self%id_mFN,   + ImFN     - mDgN_mC - mortmFN)
      _ADD_SOURCE_(self%id_mFP,   + ImFP     - mDgP_mC - mortmFP)

      total_steal = sum(self%steal)
      if (total_steal == 0) then
         _ADD_SOURCE_(self%id_mFChl, + ImFChl - mDgChl_mC - mortmFChl)
      else
         _ADD_SOURCE_(self%id_mFChl, + ImFChl - Chl_not_used - Chl_decay - mortmFChl)
      end if

      _ADD_SOURCE_(self%id_NH4,   (- mNH4up   + mregN)/14._rk)
      _ADD_SOURCE_(self%id_NO3,   - mNO3up/14._rk)
      _ADD_SOURCE_(self%id_DIP,   (- mIPu_mC  + mregP)/31._rk)
      _ADD_SOURCE_(self%id_DIC,   (- mICup    - mDOCout  + mICout)/12._rk)
      _ADD_SOURCE_(self%id_DOC,   + mDOCout/12._rk)
      _ADD_SOURCE_(self%id_sDOMC, + (mVOCL + mortC)/12._rk*self%diss_frac)
      _ADD_SOURCE_(self%id_sDOMN, + (mVONL + mortN)/14._rk*self%diss_frac)
      _ADD_SOURCE_(self%id_sDOMP, + (mVOPL + mortP)/31._rk*self%diss_frac)
      _ADD_SOURCE_(self%id_POC,   + (mVOCL + mortC)/12._rk*(1._rk - self%diss_frac))
      _ADD_SOURCE_(self%id_PON,   + (mVONL + mortN)/14._rk*(1._rk - self%diss_frac))
      _ADD_SOURCE_(self%id_POP,   + (mVOPL + mortP)/31._rk*(1._rk - self%diss_frac))

      ! Send all consumed silicate (and Si-mortality) to particulate organic matter pool.
      if (self%SSi) then 
         _ADD_SOURCE_(self%id_RPs, (sum(sI*preysP) + mortS)/28._rk)
      else
         _ADD_SOURCE_(self%id_RPs, sum(sI*preysP)/28._rk)
      end if

      _SET_DIAGNOSTIC_(self%id_mNC,     mNC)
      _SET_DIAGNOSTIC_(self%id_mPC,     mPC)
      _SET_DIAGNOSTIC_(self%id_mChlC,   mChlC)
      _SET_DIAGNOSTIC_(self%id_mFC,     mFC)
      _SET_DIAGNOSTIC_(self%id_mFChlC,  mFChlC)
      _SET_DIAGNOSTIC_(self%id_mFNC,    mFNC)
      _SET_DIAGNOSTIC_(self%id_mFPC,    mFPC)
      _SET_DIAGNOSTIC_(self%id_UmpN,    UmpN)
      _SET_DIAGNOSTIC_(self%id_mAVP,    mAVP)
      _SET_DIAGNOSTIC_(self%id_mNVP,    mNVP)
      _SET_DIAGNOSTIC_(self%id_mAV,     mAV)
      _SET_DIAGNOSTIC_(self%id_mNV,     mNV)
      _SET_DIAGNOSTIC_(self%id_mfrat,   mfrat)
      _SET_DIAGNOSTIC_(self%id_mNCu,    mNCu)
      _SET_DIAGNOSTIC_(self%id_mPCu,    mPCu)
      _SET_DIAGNOSTIC_(self%id_mNPCu,   mNPCu)
      _SET_DIAGNOSTIC_(self%id_mfrNC,   mfrNC)
      _SET_DIAGNOSTIC_(self%id_mupN,    mupN)
      _SET_DIAGNOSTIC_(self%id_mupN_mC, mupN_mC)
      _SET_DIAGNOSTIC_(self%id_mNH4up,  mNH4up)
      _SET_DIAGNOSTIC_(self%id_mNO3up,  mNO3up)
      _SET_DIAGNOSTIC_(self%id_mfrPC,   mfrPC)
      _SET_DIAGNOSTIC_(self%id_mIPu,    mIPu)
      _SET_DIAGNOSTIC_(self%id_mIPu_mC, mIPu_mC)
      _SET_DIAGNOSTIC_(self%id_PSrel,   PSrel)
      _SET_DIAGNOSTIC_(self%id_mPS,     mPS)
      _SET_DIAGNOSTIC_(self%id_mCas,    mCas)
      _SET_DIAGNOSTIC_(self%id_Rphot,   Rphot)
      _SET_DIAGNOSTIC_(self%id_Rhet,    Rhet)
      _SET_DIAGNOSTIC_(self%id_mCu,     mCu)
      _SET_DIAGNOSTIC_(self%id_mDgC,    mDgC)
      _SET_DIAGNOSTIC_(self%id_FPS,     FPS)
      _SET_DIAGNOSTIC_(self%id_mPStot,  mPStot)    
      _SET_DIAGNOSTIC_(self%id_ImF,     ImF)
      _SET_DIAGNOSTIC_(self%id_kpmax,   kpmax)
      _SET_DIAGNOSTIC_(self%id_FCmax,   FCmax)
      _SET_DIAGNOSTIC_(self%id_Cfix,    Cfix)
      _SET_DIAGNOSTIC_(self%id_Cresp,   Cresp)
      _SET_DIAGNOSTIC_(self%id_ImFChl,  ImFChl)
      _SET_DIAGNOSTIC_(self%id_ImFk,    ImFk)
      _SET_DIAGNOSTIC_(self%id_mDgChl_mC,    mDgChl_mC)
      _SET_DIAGNOSTIC_(self%id_Chl_decay,    Chl_decay)
      _SET_DIAGNOSTIC_(self%id_Chl_not_used, Chl_not_used)
      _SET_DIAGNOSTIC_(self%id_photCup, photCup)   
      _SET_DIAGNOSTIC_(self%id_hetCup,  hetCup)
      _SET_DIAGNOSTIC_(self%id_photNup, photNup)   
      _SET_DIAGNOSTIC_(self%id_hetNup,  hetNup)
      _SET_DIAGNOSTIC_(self%id_photPup, photPup)   
      _SET_DIAGNOSTIC_(self%id_hetPup,  hetPup)        
      _SET_DIAGNOSTIC_(self%id_PAR_PB,  PAR) 
      _SET_DIAGNOSTIC_(self%id_mPqm,    mPqm)
      _SET_DIAGNOSTIC_(self%id_Scon,    Scon)
      _SET_DIAGNOSTIC_(self%id_Pbal,    Pbal)
      _SET_DIAGNOSTIC_(self%id_BalRes,  BalRes)    
      _SET_DIAGNOSTIC_(self%id_PbalCon, PbalCon) 
      _SET_DIAGNOSTIC_(self%id_PD,      PD)  
      _SET_DIAGNOSTIC_(self%id_Dmax,    Dmax)  
      _SET_DIAGNOSTIC_(self%id_FCrelA,  FCrelA)     
      _SET_DIAGNOSTIC_(self%id_AE_op,   AE_op)
      _SET_DIAGNOSTIC_(self%id_Urel,    Urel)  
      _SET_DIAGNOSTIC_(self%id_mRTot,   mRTot)
      _SET_DIAGNOSTIC_(self%id_et,      et)
      _SET_DIAGNOSTIC_(self%id_mDOCout,    mDOCout)   
      _SET_DIAGNOSTIC_(self%id_mVOCL,      mVOCL) 
      _SET_DIAGNOSTIC_(self%id_mVONL,      mVONL) 
      _SET_DIAGNOSTIC_(self%id_mVOPL,      mVOPL)
      _SET_DIAGNOSTIC_(self%id_mortC,      mortC)
      _SET_DIAGNOSTIC_(self%id_mortN,      mortN)
      _SET_DIAGNOSTIC_(self%id_mortP,      mortP)
      _SET_DIAGNOSTIC_(self%id_mregN,      mregN)   
      _SET_DIAGNOSTIC_(self%id_mregP,      mregP)
      _SET_DIAGNOSTIC_(self%id_Preas,      Preas)
      _SET_DIAGNOSTIC_(self%id_NH4reas,    NH4reas)
      _SET_DIAGNOSTIC_(self%id_npp,        npp)
      _SET_DIAGNOSTIC_(self%id_FCrelV,        FCrelV)
      _SET_DIAGNOSTIC_(self%id_stress,        stress)
      _SET_DIAGNOSTIC_(self%id_mortmC,        mortmC)

      do iprey=1,self%nprey
      _SET_DIAGNOSTIC_(self%id_sI(iprey),    sI(iprey))
      _SET_DIAGNOSTIC_(self%id_Cp(iprey),    Cp(iprey))
      _SET_DIAGNOSTIC_(self%id_I(iprey),     I(iprey))
      _SET_DIAGNOSTIC_(self%id_Ccon(iprey),  Ccon(iprey))
      _SET_DIAGNOSTIC_(self%id_Icap(iprey),  Icap(iprey))
      end do

      ! Leave spatial loops (if any)
      _LOOP_END_
  end subroutine do

  subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_su_mixo), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: bc, bn, bp, bsi     
      real(rk) :: mF, mFN, mFP, mFChl, mC, mN, mP, mChl, mSi

      _BOTTOM_LOOP_BEGIN_

      ! Retrieve internal state variables values
      _GET_(self%id_mF,    mF)
      _GET_(self%id_mFN,   mFN)
      _GET_(self%id_mFP,   mFP)
      _GET_(self%id_mFChl, mFChl)
      _GET_(self%id_mC,    mC)
      _GET_(self%id_mN,    mN)
      _GET_(self%id_mP,    mP)
      _GET_(self%id_mChl,  mChl)

      _GET_BOTTOM_(self%id_bc,  bc)
      _GET_BOTTOM_(self%id_bn,  bn)
      _GET_BOTTOM_(self%id_bp,  bp)

      ! Account for sinking
      _ADD_BOTTOM_FLUX_(self%id_mC,    -self%sink * mC)
      _ADD_BOTTOM_FLUX_(self%id_mF,    -self%sink * mF)
      _ADD_BOTTOM_FLUX_(self%id_mN,    -self%sink * mN)
      _ADD_BOTTOM_FLUX_(self%id_mFN,   -self%sink * mFN)
      _ADD_BOTTOM_FLUX_(self%id_mP,    -self%sink * mP)
      _ADD_BOTTOM_FLUX_(self%id_mFP,   -self%sink * mFP)
      _ADD_BOTTOM_FLUX_(self%id_mChl,  -self%sink * mChl)
      _ADD_BOTTOM_FLUX_(self%id_mFChl, -self%sink * mFChl)

      _ADD_BOTTOM_SOURCE_(self%id_bc,         self%sink * (mC + mF)/12._rk)
      _ADD_BOTTOM_SOURCE_(self%id_bn,         self%sink * (mN + mFN)/14._rk)
      _ADD_BOTTOM_SOURCE_(self%id_bp,         self%sink * (mP + mFP)/31._rk)

      if (self%SSi) then 
      _GET_(self%id_mSi,    mSi)   
      _GET_BOTTOM_(self%id_bsi, bsi)  
      _ADD_BOTTOM_FLUX_(self%id_mSi,  -self%sink * mSi)
      _ADD_BOTTOM_SOURCE_(self%id_bsi,       self%sink * mSi/28._rk)
      end if  

      _BOTTOM_LOOP_END_
  end subroutine

end module su_mixo








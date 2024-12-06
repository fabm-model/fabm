#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_nersc_ecosmo --- ECOSMO biogeochemical model
!
! !INTERFACE:
   module fabm_nersc_ecosmo_operational
!
! !DESCRIPTION:
!
! The ECOSMO model is based on Daewel & Schrum (JMS,2013)
!
! !USES:
   use fabm_types
   use fabm_expressions
   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_nersc_ecosmo_operational
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: sedy0 = 86400.0_rk
   real(rk), parameter :: mmolm3_in_mll = 44.6608009_rk
   real(rk)            :: redf(20)=0.0_rk
   real(rk)            :: BioC(45)=0.0_rk
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_nersc_ecosmo_operational
!     Variable identifiers
      type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_sil
      type (type_state_variable_id)         :: id_opa, id_det, id_dia, id_fla
      type (type_state_variable_id)         :: id_diachl, id_flachl, id_bgchl
      type (type_state_variable_id)         :: id_mesozoo, id_microzoo, id_bg,id_dom, id_oxy
      type (type_state_variable_id)         :: id_dic, id_alk
      type (type_bottom_state_variable_id)  :: id_sed1, id_sed2, id_sed3
      type (type_dependency_id)             :: id_temp, id_salt, id_par
      type (type_dependency_id)             :: id_parmean
      type (type_horizontal_dependency_id)  :: id_tbs
      type (type_horizontal_dependency_id)  :: id_sfpar, id_meansfpar
      type (type_diagnostic_variable_id)    :: id_denit, id_primprod, id_secprod
      type (type_diagnostic_variable_id)    :: id_parmean_diag
      type (type_diagnostic_variable_id)    :: id_c2chl_fla, id_c2chl_dia,id_c2chl_bg
      type (type_diagnostic_variable_id)    :: id_nlim, id_plim, id_slim, id_llim
      type (type_horizontal_diagnostic_variable_id)    :: id_tbsout

!     Model parameters
      real(rk) :: BioC(45)
      real(rk) :: extdet, extdom
      real(rk) :: zpr, frr
      real(rk) :: prefZsPs
      real(rk) :: prefZsPl
      real(rk) :: prefZsBG
      real(rk) :: prefZsD
      real(rk) :: prefZlPs
      real(rk) :: prefZlPl
      real(rk) :: prefZlBG
      real(rk) :: prefZlD
      real(rk) :: prefZlZs
      real(rk) :: surface_deposition_no3
      real(rk) :: surface_deposition_nh4
      real(rk) :: surface_deposition_pho
      real(rk) :: surface_deposition_sil
      real(rk) :: nfixation_minimum_daily_par
      real(rk) :: bg_growth_minimum_daily_rad
      real(rk) :: MAXchl2nPs, MINchl2nPs 
      real(rk) :: MAXchl2nPl, MINchl2nPl 
      real(rk) :: MAXchl2nBG, MINchl2nBG 
      real(rk) :: alfaPl, alfaPs, alfaBG 

      ! ECOSMO modules
      logical  :: use_chl ! activates explicit chlorophyll-a variable
      logical  :: use_cyanos ! activates cyanobacteria, effective for low salinity conditions
      logical  :: couple_co2 ! activates PML (Blackford, 2004) carbon module
      logical  :: not_0d ! module to run the model in 0d
      logical  :: use_chl_in_PI_curve ! activated chl dependent light limitation (temporary now - will be permanent)
      logical  :: turn_on_additional_diagnostics ! activates additional diagnostics for model debugging 
      
      contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
      procedure :: get_light_extinction

   end type type_nersc_ecosmo_operational
!EOP
!-----------------------------------------------------------------------

   type (type_bulk_standard_variable), parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the ECOSMO model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
   class (type_nersc_ecosmo_operational),intent(inout),target  :: self
   integer,                intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
!  Caglar Yumruktepe:
!     Added fabm.yaml support: parameters from yaml file are copied to
!                           BioC array. Eventually, BioC array will be removed
!                           from the model where parameter names from the yaml
!                           file will be used.
!     Added dynamic chlorophyll-a from Geider etal., 1997
!     Added community dependent particle sinking rates
!     Added chlorophyll-a dependent light-limitation   
!
! !LOCAL VARIABLES:
! Everything else taken from yaml file
!
   integer :: i
   ! set Redfield ratios:
   redf(1) = 6.625_rk      !C_N
   redf(2) = 106.0_rk      !C_P
   redf(3) = 6.625_rk      !C_SiO
   redf(4) = 16.0_rk       !N_P
   redf(5) = 1.0_rk        !N_SiO
   redf(6) = 12.01_rk      !C_Cmg
   redf(7) = 44.6608009_rk !O2mm_ml
   redf(8) = 14.007_rk     !N_Nmg
   redf(9) = 30.97_rk      !P_Pmg
   redf(10) = 28.09_rk     !Si_Simg
   do i=1,10
     redf(i+10) = 1._rk/redf(i)
   end do


!EOP
!-----------------------------------------------------------------------
!BOC

   call self%get_parameter(self%zpr, 'zpr', '1/day', 'zpr_long_name_needed', default=0.001_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter(self%frr, 'frr', '-', 'fraction of dissolved from det.', default=0.4_rk)
   call self%get_parameter(self%nfixation_minimum_daily_par, 'nfixation_minimum_daily_par', 'nfixation minimum daily par', default=40.0_rk)
   call self%get_parameter(self%bg_growth_minimum_daily_rad, 'bg_growth_minimum_daily_rad', 'bg growth minimum daily rad', default=120.0_rk)
   ! set surface fluxes in [mgC/m2/s]
   call self%get_parameter( self%surface_deposition_no3, 'surface_deposition_no3', 'mmolN/m**2 d', 'surface deposition no3', default=0.0_rk, scale_factor=redf(1)*redf(6)/sedy0 )
   call self%get_parameter( self%surface_deposition_nh4, 'surface_deposition_nh4', 'mmolN/m**2 d', 'surface deposition nh4', default=0.0_rk, scale_factor=redf(1)*redf(6)/sedy0 )
   call self%get_parameter( self%surface_deposition_pho, 'surface_deposition_pho', 'mmolP/m**2 d', 'surface deposition pho', default=0.0_rk, scale_factor=redf(2)*redf(6)/sedy0 )
   call self%get_parameter( self%surface_deposition_sil, 'surface_deposition_sil', 'mmolSi/m**2 d', 'surface deposition sil', default=0.0_rk, scale_factor=redf(3)*redf(6)/sedy0 )
   !  change units 1/day to 1/sec and mmolN,P,Si to mmolC
   call self%get_parameter( self%BioC(1) , 'muPl',        '1/day',      'max growth rate for Pl',          default=1.30_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(2) , 'muPs',        '1/day',      'max growth rate for Ps',          default=1.10_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(3) , 'aa',          'm**2/W',     'photosynthesis ef-cy',            default=0.04_rk)
   call self%get_parameter( self%BioC(4) , 'EXw',         '1/m',        'light extinction',                default=0.041_rk)
   if (self%use_chl) then
      call self%get_parameter( self%BioC(5) , 'Exphy',       'm**2/mgCHL', 'phyto self-shading',              default=0.04_rk )
   else
      call self%get_parameter( self%BioC(5) , 'Exphy',       'm**2/mmolN', 'phyto self-shading',              default=0.04_rk, scale_factor=1.0_rk/(redf(1)*redf(6)) )
   end if
   call self%get_parameter( self%extdet , 'Exdet',       'm**2/molC', 'detritus self-shading',         default=0.0_rk, scale_factor=1.0_rk/1000.0_rk )
   call self%get_parameter( self%extdom , 'Exdom',       'm**2/molC', 'dom self-shading',              default=0.0_rk, scale_factor=1.0_rk/1000.0_rk )
   call self%get_parameter( self%BioC(6) , 'rNH4',        'mmolN/m**3', 'NH4 half saturation',             default=0.20_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%BioC(7) , 'rNO3',        'mmolN/m**3', 'NO3 half saturation',             default=0.50_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%BioC(8) , 'psi',         'm**3/mmolN', 'NH4 inhibition',                  default=3.0_rk,   scale_factor=1.0_rk/(redf(1)*redf(6)) )
   call self%get_parameter( self%BioC(9) , 'mPl',         '1/day',      'Pl mortality rate',               default=0.04_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(10), 'mPs',         '1/day',      'Ps mortality rate',               default=0.08_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(11), 'GrZlP',       '1/day',      'Grazing rate Zl on Phyto',        default=0.80_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(12), 'GrZsP',       '1/day',      'Grazing rate Zs on Phyto',        default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(13), 'GrZlZ',       '1/day',      'Grazing rate Zl on Zs',           default=0.50_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(14), 'Rg',          'mmolN/m**3', 'Zs, Zl half saturation',          default=0.50_rk,  scale_factor=redf(1)*redf(6))
   call self%get_parameter( self%BioC(15), 'mZl',         '1/day',      'Zl mortality rate',               default=0.10_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(16), 'mZs',         '1/day',      'Zs mortality rate',               default=0.20_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(17), 'excZl',       '1/day',      'Zl excretion rate',               default=0.06_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(18), 'excZs',       '1/day',      'Zs excretion rate',               default=0.08_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(19), 'gammaZlp',    '1',          'Zl assim. eff. on plankton',      default=0.75_rk)
   call self%get_parameter( self%BioC(20), 'gammaZsp',    '1',          'Zs assim. eff. on plankton',      default=0.75_rk)
   call self%get_parameter( self%BioC(21), 'gammaZd',     '1',          'Zl & Zs assim. eff. on det',      default=0.75_rk)
   call self%get_parameter( self%BioC(22), 'reminD',      '1/day',      'Detritus remin. rate',            default=0.003_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(23), 'sinkDet',     'm/day',      'Detritus sinking rate',           default=5.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(24), 'Wa',          'm/day',      '???',                             default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(25),  'rPO4',       'mmolP/m**3', 'PO4 half saturation',             default=0.05_rk,  scale_factor=redf(2)*redf(6))
   call self%get_parameter( self%BioC(26),  'rSi',        'mmolSi/m**3','SiO2 half saturation',            default=0.50_rk,  scale_factor=redf(3)*redf(6))
   call self%get_parameter( self%BioC(27),  'regenSi',    '1/day',      'Si regeneration rate',            default=0.015_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(28),  'muBG',       '1/day',      'max growth rate for BG',          default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(29),  'TctrlBG',    '1/degC',     'BG T control beta',               default=1.00_rk)
   call self%get_parameter( self%BioC(30),  'TrefBG',     'degC',       'BG reference temperature',        default=0.00_rk)
   call self%get_parameter( self%BioC(31),  'GrBG',       '1/day',      'BG max grazing rate',             default=0.30_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(32),  'mBG',        '1/day',      'BG mortality rate',               default=0.08_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(33),  'upliftBG',   'm/day',      'BG uplifting rate',               default=0.10_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(34),  'crBotStr',   'N/m**2',     'critic. bot. stress for resusp.', default=0.007_rk)
   call self%get_parameter( self%BioC(35),  'resuspRt',   '1/day',      'resuspension rate',               default=25.0_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(36),  'sedimRt',    'm/day',      'sedimentation rate',              default=3.5_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(37),  'burialRt',   '1/day',      'burial rate',                     default=1E-5_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(38),  'reminSED',   '1/day',      'sediment remineralization rate',  default=0.001_rk, scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(39),  'TctrlDenit', '1/degC',     'temp. control denitrification',   default=0.15_rk)
   call self%get_parameter( self%BioC(40),  'RelSEDp1',   'units??',    'P sedim. rel. p1',                default=0.15_rk)
   call self%get_parameter( self%BioC(41),  'RelSEDp2',   'units??',    'P sedim. rel. p2',                default=0.10_rk)
   call self%get_parameter( self%BioC(42),  'reminSEDsi', '1/day',      'sed. remineralization rate Si',   default=0.0002_rk,scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(43),  'sinkOPAL',   'm/day',      'OPAL sinking rate',               default=5.0_rk,   scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(44),  'sinkBG',     'm/day',      'BG sinking rate',                 default=-1.0_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%BioC(45),  'sinkDia',    'm/day',      'Diatom sinking rate',             default=0.0_rk,   scale_factor=1.0_rk/sedy0)
   !  growth fractions
   call self%get_parameter( self%prefZsPs,  'prefZsPs',   '-',          'Grazing preference Zs on Ps',     default=0.70_rk)
   call self%get_parameter( self%prefZsPl,  'prefZsPl',   '-',          'Grazing preference Zs on Pl',     default=0.25_rk)
   call self%get_parameter( self%prefZsD,   'prefZsD',    '-',          'Grazing preference Zs on Det.',   default=0.00_rk)
   call self%get_parameter( self%prefZsBG,  'prefZsBG',   '-',          'Grazing preference Zs on BG',     default=0.30_rk)
   call self%get_parameter( self%prefZlPs,  'prefZlPs',   '-',          'Grazing preference Zl on Ps',     default=0.10_rk)
   call self%get_parameter( self%prefZlPl,  'prefZlPl',   '-',          'Grazing preference Zl on Pl',     default=0.85_rk)
   call self%get_parameter( self%prefZlZs,  'prefZlZs',   '-',          'Grazing preference Zl on Zs',     default=0.15_rk)
   call self%get_parameter( self%prefZlD,   'prefZlD',    '-',          'Grazing preference Zl on Det.',   default=0.00_rk)
   call self%get_parameter( self%prefZlBG,  'prefZlBG',   '-',          'Grazing preference Zl on BG',     default=0.30_rk)
   ! chlorophyll-a constants
   call self%get_parameter( self%MINchl2nPs, 'MINchl2nPs', 'mgChl/mmolN', 'minimum Chl to N ratio Ps', default=0.50_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MAXchl2nPs, 'MAXchl2nPs', 'mgChl/mmolN', 'maximum Chl to N ratio Ps', default=3.83_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MINchl2nPl, 'MINchl2nPl', 'mgChl/mmolN', 'minimum Chl to N ratio Pl', default=0.50_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MAXchl2nPl, 'MAXchl2nPl', 'mgChl/mmolN', 'maximum Chl to N ratio Pl', default=2.94_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MINchl2nBG, 'MINchl2nBG', 'mgChl/mmolN', 'minimum Chl to N ratio BG', default=0.50_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%MAXchl2nBG, 'MAXchl2nBG', 'mgChl/mmolN', 'maximum Chl to N ratio BG', default=3.83_rk, scale_factor=redf(11)*redf(16))
   call self%get_parameter( self%alfaPs,     'alfaPs', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve Ps', default=0.0393_rk, scale_factor=redf(1)*redf(6) )
   call self%get_parameter( self%alfaPl,     'alfaPl', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve Pl', default=0.0531_rk, scale_factor=redf(1)*redf(6) )
   call self%get_parameter( self%alfaBG,     'alfaBG', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve BG', default=0.0393_rk, scale_factor=redf(1)*redf(6) )
   ! add switches
   call self%get_parameter( self%use_cyanos,     'use_cyanos', '', 'switch cyanobacteria', default=.true.)
   call self%get_parameter( self%couple_co2,     'couple_co2', '', 'switch coupling to carbonate module', default=.false.)
   call self%get_parameter( self%use_chl,     'use_chl', '', 'switch chlorophyll/c dynamics', default=.true.)
   call self%get_parameter( self%not_0d,     'not_0d', '', 'do not run the model in a 0D box', default=.true.)
   call self%get_parameter( self%use_chl_in_PI_curve, 'use_chl_in_PI_curve','','activated chl dependent light limitation',default=.false.)
   call self%get_parameter( self%turn_on_additional_diagnostics, 'turn_on_additional_diagnostics','','activates additional diagnostics for model debugging',default=.false.)
   ! Register state variables
   call self%register_state_variable( self%id_no3,      'no3',     'mgC/m3',    'nitrate',                   minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=5.0_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_nh4,      'nh4',     'mgC/m3',    'ammonium',                  minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=0.1_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_pho,      'pho',     'mgC/m3',    'phosphate',                 minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=0.3_rk*redf(2)*redf(6)  )
   call self%register_state_variable( self%id_sil,      'sil',     'mgC/m3',    'silicate',                  minimum=0.0_rk,        vertical_movement=0.0_rk,  &
                                      initial_value=5.0_rk*redf(3)*redf(6)  )
   call self%register_state_variable( self%id_oxy,      'oxy',     'mmolO2/m3', 'oxygen',                    minimum=0.0_rk,   vertical_movement=0.0_rk,  &
                                      initial_value=85.0_rk  )
   call self%register_state_variable( self%id_fla,      'fla',     'mgC/m3',    'small phytoplankton',       minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-4_rk*redf(1)*redf(6))
   call self%register_state_variable( self%id_dia,      'dia',     'mgC/m3',    'large phytoplankton',       minimum=1.0e-7_rk,     vertical_movement=-self%BioC(45) , &
                                      initial_value=1e-4_rk*redf(1)*redf(6) )
   if (self%use_cyanos) then
     call self%register_state_variable( self%id_bg,       'bg',      'mgC/m3',    'cyanobacteria',             minimum=1.0e-14_rk,     vertical_movement=-self%BioC(44) , &
                                      initial_value=1e-4_rk*redf(1)*redf(6) )
     if (self%use_chl) then
       call self%register_state_variable( self%id_bgchl,    'bgchl',   'mgChl/m3',  'cyanobacteria chl-a',       minimum=1.0e-14_rk/20., vertical_movement=-self%BioC(44) , &
                                      initial_value=1e-4_rk*redf(1)*redf(6)/27.)
       call self%add_to_aggregate_variable(total_chlorophyll, self%id_bgchl)
     else
       call self%add_to_aggregate_variable(total_chlorophyll, self%id_bg, scale_factor=1.0_rk/60.0_rk)
     end if
   end if
   if (self%use_chl) then
     call self%register_state_variable( self%id_diachl,   'diachl',  'mgChl/m3',  'large phytoplankton chl-a', minimum=1.0e-7_rk/27., vertical_movement=-self%BioC(45) , &
                                      initial_value=1e-4_rk*redf(1)*redf(6)/27.)
     call self%register_state_variable( self%id_flachl,   'flachl',  'mgChl/m3',  'small phytoplankton chl-a', minimum=1.0e-7_rk/20., vertical_movement=0.0_rk, &
                                      initial_value=1e-4_rk*redf(1)*redf(6)/20.)
     call self%add_to_aggregate_variable(total_chlorophyll, self%id_diachl)
     call self%add_to_aggregate_variable(total_chlorophyll, self%id_flachl)
   else
    call self%add_to_aggregate_variable(total_chlorophyll, self%id_dia, scale_factor=1.0_rk/60.0_rk)
    call self%add_to_aggregate_variable(total_chlorophyll, self%id_fla, scale_factor=1.0_rk/60.0_rk)
   end if
   call self%register_state_variable( self%id_microzoo, 'microzoo','mgC/m3',    'microzooplankton',          minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-6_rk*redf(1)*redf(6) )
   call self%register_state_variable( self%id_mesozoo,  'mesozoo', 'mgC/m3',    'mesozooplankton',           minimum=1.0e-7_rk,     vertical_movement=0.0_rk, &
                                      initial_value=1e-6_rk*redf(1)*redf(6) )
   call self%register_state_variable( self%id_det,      'det',     'mgC/m3',    'detritus',                  minimum=0.0_rk,   vertical_movement=-self%BioC(23),      &
                                      initial_value=2.0_rk*redf(1)*redf(6)  )
   call self%register_state_variable( self%id_opa,      'opa',     'mgC/m3',    'opal',                      minimum=0.0_rk, vertical_movement=-self%BioC(43), &
                                      initial_value=2.0_rk*redf(3)*redf(6) )
   call self%register_state_variable( self%id_dom,      'dom',     'mgC/m3',    'labile dissolved om',       minimum=0.0_rk , &
                                      initial_value=3.0_rk*redf(1)*redf(6)   )
   call self%register_state_variable( self%id_sed1,     'sed1',    'mgC/m2',    'sediment detritus',         minimum=0.0_rk , &
                                      initial_value=20.0_rk*redf(1)*redf(6)*redf(18) )
   call self%register_state_variable( self%id_sed2,     'sed2',    'mgC/m2',    'sediment opal',             minimum=0.0_rk , &
                                      initial_value=20.0_rk*redf(3)*redf(6)*redf(20) )
   call self%register_state_variable( self%id_sed3,     'sed3',    'mgC/m2',    'sediment adsorbed pho.',    minimum=0.0_rk , &
                                      initial_value=2.0_rk*redf(2)*redf(6)*redf(19) )
   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_primprod,'primprod','mgC/m**3/s', &
         'primary production rate', output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_secprod,'secprod','mgC/m**3/s', &
         'secondary production rate', output=output_time_step_averaged)

   if (self%turn_on_additional_diagnostics) then  
         call self%register_diagnostic_variable(self%id_parmean_diag,'parmean','W/m**2', &
         'daily-mean photosynthetically active radiation', output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_nlim,'Nlim','-', &
         'N-limitation on primary production', output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_plim,'Plim','-', &
         'P-limitation on primary production', output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_slim,'Slim','-', &
         'Si-limitation on primary production', output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_llim,'Llim','-', &
         'Light-limitation on primary production', output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_denit,'denit','mmolN/m**3/s', &
         'denitrification rate', output=output_time_step_averaged)
         call self%register_diagnostic_variable(self%id_tbsout,'botstrss','fill_later', &
         'total bottom stress', source=source_do_bottom)
         if (self%use_chl) then
          call self%register_diagnostic_variable(self%id_c2chl_fla,'c2chl_fla','mgC/mgCHL', &
              'daily-mean C to CHL ratio for flagellates', output=output_time_step_averaged)
          call self%register_diagnostic_variable(self%id_c2chl_dia,'c2chl_dia','mgC/mgCHL', &
              'daily-mean C to CHL ratio for diatoms', output=output_time_step_averaged)
              if (self%use_cyanos) then
                  call self%register_diagnostic_variable(self%id_c2chl_bg,'c2chl_bg','mgC/mgCHL', &
                  'daily-mean C to CHL ratio for cyanobacteria', output=output_time_step_averaged)
              end if
         end if
   end if

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   if (self%not_0d) then
      call self%register_dependency(self%id_tbs,standard_variables%bottom_stress)
   end if 
   call self%register_dependency(self%id_sfpar,standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   ! use temporal mean of light for the last 24 hours
   call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk))
   call self%register_dependency(self%id_meansfpar,temporal_mean(self%id_sfpar,period=86400._rk,resolution=3600._rk))

   if (self%couple_co2) then
     call self%register_state_dependency(self%id_dic, 'dic_target','mmol m-3','dic budget')
     call self%register_state_dependency(self%id_alk, 'alk_target','mmol m-3','alkalinity budget')
   end if

!   call self%register_dependency(self%id_h,       'icethickness', 'm',    'ice thickness')
!   call self%register_dependency(self%id_hs,      'snowthickness','m',    'snow thickness')
   

   return

end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of ECOSMO model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_nersc_ecosmo_operational),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   real(rk) :: no3,nh4,pho,sil,t_sil,oxy,fla,bg,dia
   real(rk) :: flachl,diachl,bgchl,chl2c_fla,chl2c_dia,chl2c_bg
   real(rk) :: microzoo,mesozoo,opa,det,dom
   real(rk) :: temp,salt,par
   real(rk) :: frem, fremDOM, blight
   real(rk) :: blightDIA,blightFLA,blightCYA ! community specific P-I curves
   real(rk) :: Ts,Tl,Tbg
   real(rk) :: Prod,Ps_prod,Pl_prod,Bg_prod
   real(rk) :: Fs,Fl,ZlonPs,ZlonPl,ZsonD,ZlonD,ZlonBg,ZsonBg,ZsonPs,ZsonPl,ZlonZs
   real(rk) :: up_no3,up_nh4,up_n,up_pho,up_sil
   real(rk) :: bioom1,bioom2,bioom3,bioom4,bioom5,bioom6,bioom7,bioom8,Onitr
   real(rk) :: rhs,dxxdet
   real(rk) :: Zl_prod, Zs_prod
   real(rk) :: mean_par, mean_surface_par, Bg_fix
   real(rk) :: fla_loss=1.0_rk
   real(rk) :: dia_loss=1.0_rk
   real(rk) :: bg_loss=1.0_rk
   real(rk) :: mic_loss=1.0_rk
   real(rk) :: mes_loss=1.0_rk
   real(rk) :: tbs
   real(rk) :: rhs_oxy,rhs_amm,rhs_nit
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_(self%id_par,par)
   _GET_(self%id_no3,no3)
   _GET_(self%id_nh4,nh4)
   _GET_(self%id_pho,pho)
   _GET_(self%id_sil,sil)
   _GET_(self%id_dia,dia)
   _GET_(self%id_fla,fla)
   if (self%use_cyanos) then
     _GET_(self%id_bg,bg)
   else
     bg=0.0_rk
   end if
   if (self%use_chl) then
     _GET_(self%id_diachl,diachl)
     _GET_(self%id_flachl,flachl)
     if (self%use_cyanos) then
       _GET_(self%id_bgchl,bgchl)
     end if
   end if
   _GET_(self%id_microzoo,microzoo)
   _GET_(self%id_mesozoo,mesozoo)
   _GET_(self%id_det,det)
   _GET_(self%id_dom,dom)
   _GET_(self%id_opa,opa)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_parmean,mean_par)
   _GET_HORIZONTAL_(self%id_meansfpar,mean_surface_par)
   _GET_HORIZONTAL_(self%id_tbs,tbs)

   ! CAGLAR
   ! checks - whether the biomass of plankton is below a predefined threshold,
   !          where below the threshold, loss terms are removed from the RHS of
   !          the equations. The idea is to keep plankton safe from extinction.
   ! loss terms are multiplied by the constants below, which can only be set
   ! by the model to 0 or 1.

   fla_loss = max(sign(-1.0_rk,fla-0.1_rk),0.0_rk)       ! flagellates
   dia_loss = max(sign(-1.0_rk,dia-0.1_rk),0.0_rk)       ! diatoms
   bg_loss  = max(sign(-1.0_rk,bg-0.1_rk),0.0_rk)        ! cyanobacteria
   mic_loss = max(sign(-1.0_rk,microzoo-0.01_rk),0.0_rk) !microzooplankton
   mes_loss = max(sign(-1.0_rk,mesozoo-0.01_rk),0.0_rk) ! mesozooplankton

   ! remineralisation rate
   frem = self%BioC(22) * (1._rk+20._rk*(temp**2/(13._rk**2+temp**2)))
   fremDOM = 10.0_rk * frem

   ! nutrient limitation factors
   ! k denotes half-saturation values
   up_nh4 = nh4/(self%BioC(6)+nh4)
   up_no3 = no3/(self%BioC(7)+no3)*exp(-self%BioC(8)*nh4)
   up_n = up_nh4+up_no3
   up_pho = pho/(self%BioC(25)+pho)
   t_sil = max(sil-80._rk,0.0_rk)
   up_sil = t_sil/(self%BioC(26)+t_sil)

   ! temperature dependence
   Ts = 1.0_rk
   Tl = 1.0_rk
   if ((salt<=10.0) .and. (mean_surface_par > self%bg_growth_minimum_daily_rad)) then
     Tbg = 1.0_rk/(1.0_rk + exp(self%BioC(29)*(self%BioC(30)-temp)))
   else
     Tbg = 0.0_rk
   end if

   ! light limitation - old formulation - default for now - CAGLAR:Sep2022
   blight = max(tanh(self%BioC(3)*par),0.0_rk)
   Ps_prod = Ts * min(blight, up_n, up_pho)
   Pl_prod = Tl * min(blight, up_n, up_pho, up_sil)
   ! production and nutrient uptake
   if (self%use_chl_in_PI_curve) then 
      blightFLA  = max( ((flachl/fla)*self%alfaPs*par) / sqrt((self%BioC(2)*sedy0)**2+(flachl/fla)**2*self%alfaPs**2*(par**2)) ,0.0_rk)
      blightDIA  = max( ((diachl/dia)*self%alfaPl*par) / sqrt((self%BioC(1)*sedy0)**2+(diachl/dia)**2*self%alfaPl**2*(par**2)) ,0.0_rk)
      Ps_prod = Ts * min(blightFLA, up_n, up_pho)
      Pl_prod = Tl * min(blightDIA, up_n, up_pho, up_sil)
   end if
   
   Prod = self%BioC(1)*Pl_prod*dia + & ! diatoms production
          self%BioC(2)*Ps_prod*fla     ! flagellates production

   if (self%use_cyanos) then
     Bg_prod = Tbg * min(blight, up_n, up_pho)
     Prod = Prod + self%BioC(28)*Bg_prod*bg ! cyanobacteria production
   end if

   if (self%use_cyanos .and. mean_par > self%nfixation_minimum_daily_par) then
     Bg_fix = Tbg * min(blight, up_pho) - Bg_prod
   else
     Bg_fix = 0.0_rk
   end if

   if (self%use_chl) then   
     ! chlorophyll-a to C change
     chl2c_fla = self%MAXchl2nPs * max(0.01_rk,Ps_prod) * self%BioC(2) * sedy0 * fla / &
               (self%alfaPs * par * flachl)
     chl2c_dia = self%MAXchl2nPl * max(0.01_rk,Pl_prod) * self%BioC(1) * sedy0 * dia / &
               (self%alfaPl * par * diachl)

            chl2c_fla = max(self%MINchl2nPs,chl2c_fla)
            chl2c_fla = min(self%MAXchl2nPs,chl2c_fla)
            chl2c_dia = max(self%MINchl2nPl,chl2c_dia)
            chl2c_dia = min(self%MAXchl2nPl,chl2c_dia)
     if (self%use_cyanos) then
       chl2c_bg = self%MAXchl2nBG * max(0.1_rk,Bg_prod) * self%BioC(28) * sedy0 * bg / &
               (self%alfaBG * par * bgchl)
       chl2c_bg  = max(self%MINchl2nBG,chl2c_bg)
       chl2c_bg  = min(self%MAXchl2nBG,chl2c_bg)
     end if
   end if


! EXPERIMENTING WITH ZOOPLANKTON GRAZING PRESSURES
   ! grazing
!   Fs = self%prefZsPs*fla*(fla/(fla+12.)) + self%prefZsPl*dia*(dia/(dia+12.)) + self%prefZsD*det*(det/(det+12.)) + self%prefZsBG*bg*(bg/(bg+12.))
!   Fl = self%prefZlPs*fla*(fla/(fla+12.)) + self%prefZlPl*dia*(dia/(dia+12.)) + self%prefZlZs*microzoo*(microzoo/(microzoo+12.)) + &
!         self%prefZlD*det*(det/(det+12.)) + self%prefZlBG*bg*(bg/(bg+12.))
!   if (self%use_cyanos) then
!    Fs = Fs + self%prefZsBG*bg
!    Fl = Fl + self%prefZlBg*bg
!   end if
!
!   ZsonPs = fla_loss * self%BioC(12) * self%prefZsPs * fla*(fla/(fla+12.))/(self%BioC(14) + Fs)
!   ZsonPl = dia_loss * self%BioC(12) * self%prefZsPl * dia*(dia/(dia+12.))/(self%BioC(14) + Fs)
!   ZsonD  =            self%BioC(12) * self%prefZsD * det*(det/(det+12.))/(self%BioC(14) + Fs)
!
!   ZlonPs = fla_loss * self%BioC(11) * self%prefZlPs * fla*(fla/(fla+12.))/(self%BioC(14) + Fl)
!   ZlonPl = dia_loss * self%BioC(11) * self%prefZlPl * dia*(dia/(dia+12.))/(self%BioC(14) + Fl)
!   ZlonD =             self%BioC(11) * self%prefZlD * det*(det/(det+12.))/(self%BioC(14) + Fl)
!   ZlonZs = mic_loss * self%BioC(13) * self%prefZlZs * microzoo*(microzoo/(microzoo+12.))/(self%BioC(14) + Fl)
!   if (self%use_cyanos) then
!     ZsonBg = bg_loss  * self%BioC(31) * self%prefZsBG * bg/(self%BioC(14) + Fs)
!     ZlonBg = bg_loss  * self%BioC(31) * self%prefZlBG * bg/(self%BioC(14) + Fl)
!   else
!     ZsonBg=0.0_rk
!     ZlonBg=0.0_rk
!   end if
! EXPERIMENTAL

   Fs = self%prefZsPs*fla + self%prefZsPl*dia + self%prefZsD*det + self%prefZsBG*bg
   Fl = self%prefZlPs*fla + self%prefZlPl*dia + self%prefZlZs*microzoo + &
         self%prefZlD*det + self%prefZlBG*bg
   if (self%use_cyanos) then
    Fs = Fs + self%prefZsBG*bg
    Fl = Fl + self%prefZlBg*bg
   end if

   ZsonPs = fla_loss * self%BioC(12) * self%prefZsPs * fla/(self%BioC(14) + Fs)
   ZsonPl = dia_loss * self%BioC(12) * self%prefZsPl * dia/(self%BioC(14) + Fs)
   ZsonD  =            self%BioC(12) * self%prefZsD * det/(self%BioC(14) + Fs)

   ZlonPs = fla_loss * self%BioC(11) * self%prefZlPs * fla/(self%BioC(14) + Fl)
   ZlonPl = dia_loss * self%BioC(11) * self%prefZlPl * dia/(self%BioC(14) + Fl)
   ZlonD =             self%BioC(11) * self%prefZlD * det/(self%BioC(14) + Fl)
   ZlonZs = mic_loss * self%BioC(13) * self%prefZlZs * microzoo/(self%BioC(14) + Fl)

   ! EXPERIMENT
   !ZsonPs = fla_loss * self%BioC(12) * self%prefZsPs * fla**2/(self%BioC(14)**2 + Fs**2)
   !ZsonPl = dia_loss * self%BioC(12) * self%prefZsPl * dia**2/(self%BioC(14)**2 + Fs**2)
   !ZsonD  =            self%BioC(12) * self%prefZsD * det**2/(self%BioC(14)**2 + Fs**2)

   !ZlonPs = fla_loss * self%BioC(11) * self%prefZlPs * fla**2/(self%BioC(14)**2 + Fl**2)
   !ZlonPl = dia_loss * self%BioC(11) * self%prefZlPl * dia**2/(self%BioC(14)**2 + Fl**2)
   !ZlonD =             self%BioC(11) * self%prefZlD * det**2/(self%BioC(14)**2 + Fl**2)
   !ZlonZs = mic_loss * self%BioC(13) * self%prefZlZs * microzoo**2/(self%BioC(14)**2 + Fl**2)
   ! EXPERIMENT

   if (self%use_cyanos) then
     ZsonBg = bg_loss  * self%BioC(31) * self%prefZsBG * bg/(self%BioC(14) + Fs)
     ZlonBg = bg_loss  * self%BioC(31) * self%prefZlBG * bg/(self%BioC(14) + Fl)
   else
     ZsonBg=0.0_rk
     ZlonBg=0.0_rk
   end if


   ! nitrification
   Onitr = 0.01_rk * redf(7) !according to Neumann  (Onitr in mlO2/l see also Stigebrand and Wulff)
   bioom1 = 0.0_rk
   bioom2 = 0.0_rk
   bioom3 = 0.0_rk
   bioom4 = 0.0_rk
   bioom5 = 0.0_rk
   bioom6 = 0.0_rk
   bioom7 = 0.0_rk
   bioom8 = 0.0_rk
   if (oxy > 0) then
     bioom1 = 0.1_rk/secs_pr_day * exp(temp*0.11_rk) * oxy/(Onitr+oxy)
     bioom2 = bioom1
     bioom6 = 1.0_rk
   else
     if (no3>0) then
       bioom5 = 5.0_rk
       bioom8 = 1.0_rk
     else
       bioom7 = 1.0_rk
     end if
   end if

! reaction rates

   _ADD_SOURCE_(self%id_fla, (self%BioC(2)*Ps_prod - self%BioC(10)*fla_loss)*fla - ZsonPs*microzoo - ZlonPs*mesozoo)
   _ADD_SOURCE_(self%id_dia, (self%BioC(1)*Pl_prod - self%BioC(9)*dia_loss)*dia - ZsonPl*microzoo - ZlonPl*mesozoo)
   if (self%use_cyanos) then
     _ADD_SOURCE_(self%id_bg,  (self%BioC(28)*(Bg_prod + Bg_fix) - self%BioC(32)*bg_loss)*bg - ZsonBg*microzoo - ZlonBg*mesozoo)
   end if

  ! for chlorophyll-a
   if (self%use_chl) then
     rhs = self%BioC(2)*Ps_prod*chl2c_fla*fla - ( (self%BioC(10)*fla_loss*fla + ZsonPs*microzoo + ZlonPs*mesozoo)*flachl/fla )
     _ADD_SOURCE_(self%id_flachl,rhs)
     rhs = self%BioC(1)*Pl_prod*chl2c_dia*dia - ( (self%BioC(9)*dia*dia_loss + ZsonPl*microzoo + ZlonPl*mesozoo)*diachl/dia )
     _ADD_SOURCE_(self%id_diachl,rhs)
     if (self%use_cyanos) then
       rhs = self%BioC(28)*(Bg_prod + Bg_fix)*chl2c_bg*bg - ((self%BioC(32)*bg*bg_loss + ZsonBg*microzoo + ZlonBg*mesozoo)*bgchl/bg )
       _ADD_SOURCE_(self%id_bgchl,rhs)
     end if  
   end if

   ! microzooplankton

   Zs_prod = self%BioC(20)*(ZsonPs + ZsonPl + ZsonBg) + self%BioC(21)*ZsonD
   rhs = (Zs_prod - (self%BioC(16) + self%BioC(18) + self%zpr)*mic_loss) * microzoo &
         - ZlonZs * mesozoo
   _ADD_SOURCE_(self%id_microzoo, rhs)

   ! mesozooplankton
   Zl_prod = self%BioC(19)*(ZlonPs + ZlonPl + ZlonBg + ZlonZs) + self%BioC(21)*ZlonD
   rhs = (Zl_prod - (self%BioC(15) + self%BioC(17) + self%zpr)*mes_loss) * mesozoo
   _ADD_SOURCE_(self%id_mesozoo, rhs)

   ! detritus
   dxxdet = (  ((1.0_rk-self%BioC(20))*(ZsonPs + ZsonPl + ZsonBg) &
              + (1.0_rk-self%BioC(21))*ZsonD) * microzoo &
              + ((1.0_rk-self%BioC(19))*(ZlonPs + ZlonPl + ZlonBg + ZlonZs) &
              + (1.0_rk-self%BioC(21))*ZlonD) * mesozoo &
              + self%BioC(16) * microzoo * mic_loss &
              + self%BioC(15) * mesozoo * mes_loss &
              + self%BioC(10) * fla * fla_loss &
              + self%BioC(9)  * dia * dia_loss)
    if (self%use_cyanos) then
       dxxdet = dxxdet + (self%BioC(32) * bg * bg_loss )
    end if

   rhs = (1.0_rk-self%frr) * dxxdet &
         - ZsonD * microzoo &
         - ZlonD * mesozoo &
         - frem * det
   _ADD_SOURCE_(self%id_det, rhs)

   ! labile dissolved organic matter
   _ADD_SOURCE_(self%id_dom, self%frr*dxxdet - fremDOM * dom)

   ! nitrate
   rhs_nit = -(up_no3+0.5d-10)/(up_n+1.0d-10)*Prod &
         + bioom1 * nh4 &
         - bioom3 * no3 &
         - frem * det * bioom5 &
         - fremDOM * dom * bioOM5
   _ADD_SOURCE_(self%id_no3, rhs_nit)

   ! ammonium
   rhs_amm = -(up_nh4+0.5d-10)/(up_n+1.0d-10)*Prod &
         + self%BioC(18) * microzoo * mic_loss &
         + self%BioC(17) * mesozoo * mes_loss &
         + frem * det &
         + fremDOM * dom - bioom1 * nh4
   _ADD_SOURCE_(self%id_nh4, rhs_amm)

   ! phosphate

   rhs = -Prod -self%BioC(28)*bg*Bg_fix &
         + self%BioC(18) * microzoo * mic_loss &
         + self%BioC(17) * mesozoo * mes_loss &
         + frem*det + fremDOM*dom
   _ADD_SOURCE_(self%id_pho, rhs)


   ! silicate
   _ADD_SOURCE_(self%id_sil, -self%BioC(1)*Pl_prod*dia + self%BioC(27)*opa)

   ! opal
   _ADD_SOURCE_(self%id_opa, self%BioC(9)*dia*dia_loss + ZsonPl*microzoo + ZlonPl*mesozoo - self%BioC(27)*opa)

   ! oxygen
   rhs_oxy = ((6.625*up_nh4 + 8.125*up_no3+1.d-10)/(up_n+1.d-10)*Prod &
         -bioom6*6.625*(self%BioC(18)*microzoo*mic_loss &
         +self%BioC(17)*mesozoo*mes_loss) &
         -frem*det*(bioom6+bioom7)*6.625 &
         -(bioom6+bioom7)*6.625*fremDOM*dom &
         -2.0_rk*bioom1*nh4)*redf(11)*redf(16)
   _ADD_SOURCE_(self%id_oxy, rhs_oxy)

   ! Carbonate dynamics
   if (self%couple_co2) then
     rhs =  redf(16) *( self%BioC(18)*microzoo*mic_loss &
            + self%BioC(17)*mesozoo*mes_loss &
            + frem*det + fremDOM*dom - Prod) 
     _ADD_SOURCE_(self%id_dic, rhs)

     rhs = redf(16)*(rhs_amm-rhs_nit)*redf(11) - 0.5_rk*rhs_oxy*(1._rk-bioom6)
     _ADD_SOURCE_(self%id_alk, rhs)
   end if

   ! Export diagnostic variables
   
   _SET_DIAGNOSTIC_(self%id_primprod, Prod + self%BioC(28)*bg*Bg_fix )
   _SET_DIAGNOSTIC_(self%id_secprod, Zl_prod*mesozoo + Zs_prod*microzoo)

   if (self%turn_on_additional_diagnostics) then
    _SET_DIAGNOSTIC_(self%id_parmean_diag, mean_par)
    _SET_DIAGNOSTIC_(self%id_denit,(frem*det*bioom5+fremDOM*dom*bioom5)*redf(11)*redf(16))
    
    _SET_DIAGNOSTIC_(self%id_nlim, up_n)
    _SET_DIAGNOSTIC_(self%id_plim, up_pho)
    _SET_DIAGNOSTIC_(self%id_slim, up_sil)
    if (self%use_chl_in_PI_curve) then
      _SET_DIAGNOSTIC_(self%id_llim, (blightFLA + blightDIA)/2.) ! cyanobacteria is not included here
    else
      _SET_DIAGNOSTIC_(self%id_llim, blight)
    end if


    if (self%use_chl) then
       _SET_DIAGNOSTIC_(self%id_c2chl_fla, 1.0_rk/chl2c_fla)
       _SET_DIAGNOSTIC_(self%id_c2chl_dia, 1.0_rk/chl2c_dia)
       if (self%use_cyanos) then
         _SET_DIAGNOSTIC_(self%id_c2chl_bg, 1.0_rk/chl2c_bg)
       end if
    end if
   end if
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the ecosmo model
!
! !INTERFACE:

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_nersc_ecosmo_operational),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk) :: o2flux, T, tr, S, o2sat, oxy
   real(rk) :: no3flux, phoflux
   real(rk) :: pho,par,bg,blight,tbg,up_pho,prod
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,T)
   _GET_(self%id_salt,S)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_par,par)
   _GET_(self%id_pho,pho)
   if (self%use_cyanos) then
     _GET_(self%id_bg,bg)
   else
     bg=0.0_rk
   end if

   ! Oxygen saturation micromol/liter__(Benson and Krause, 1984)
   tr = 1.0_rk/(T + 273.15_rk)
   o2sat= exp(- 135.90205_rk              &
       + (1.575701d05 ) * tr               &
       - (6.642308d07 ) * tr**2            &
       + (1.243800d10) * tr**3            &
       - (8.621949d11) * tr**4            &
       - S*(0.017674_rk-10.754_rk*tr+2140.7_rk*tr**2)  )

!   o2flux = 5._rk/secs_pr_day * (o2sat - oxy)
   o2flux = 1._rk/secs_pr_day * (o2sat - oxy)

   _ADD_SURFACE_FLUX_(self%id_oxy,o2flux)

   _ADD_SURFACE_FLUX_(self%id_no3,self%surface_deposition_no3)
   _ADD_SURFACE_FLUX_(self%id_nh4,self%surface_deposition_nh4)
   _ADD_SURFACE_FLUX_(self%id_pho,self%surface_deposition_pho)
   _ADD_SURFACE_FLUX_(self%id_sil,self%surface_deposition_sil)

#if 0
   if (self%use_cyanos) then
     ! calculate cyanobacteria surface production
     if (S <= 10.0) then
       tbg = 1.0_rk/(1.0_rk + exp(self%BioC(29)*(self%BioC(30)-T)))
     else
       tbg = 0.0_rk
     end if

     blight=max(tanh(self%BioC(3)*par),0.)
     up_pho = pho/(self%BioC(25)+pho)
     prod = self%BioC(28) * bg * Tbg * min(blight, up_pho) ! cyanobacteria production

   !_ADD_SOURCE_(self%id_bg,  prod)
   !_ADD_SOURCE_(self%id_pho, -prod)
   !_ADD_SURFACE_SOURCE_(id_oxy, ) ! not included in the modular ECOSMO version
   !_ADD_SURFACE_SOURCE_(id_dic, -Prod)
   end if
#endif

   ! Leave spatial loops over the horizontal domain (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Bottom fluxes for the ecosmo model
!
! !INTERFACE:

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_nersc_ecosmo_operational),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk) :: temp, tbs, oxy, no3, det, opa, sed1, sed2, sed3
   real(rk) :: pho, Rds, Rsd, Rsa, Rsdenit, Rsa_p, yt1, yt2
   real(rk) :: rhs, flux, alk_flux
   real(rk) :: bioom1, bioom2, bioom3, bioom4, bioom5, bioom6, bioom7, bioom8
   ! add community sinking local variables
   real(rk) :: dsnk
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_det,det)
   _GET_(self%id_opa,opa)
   _GET_(self%id_no3,no3)
   _GET_HORIZONTAL_(self%id_sed1,sed1)
   _GET_HORIZONTAL_(self%id_sed2,sed2)
   _GET_HORIZONTAL_(self%id_sed3,sed3)
   _GET_HORIZONTAL_(self%id_tbs,tbs)

   bioom1 = 0.0_rk
   bioom2 = 0.0_rk
   bioom3 = 0.0_rk
   bioom4 = 0.0_rk
   bioom5 = 0.0_rk
   bioom6 = 0.0_rk
   bioom7 = 0.0_rk
   bioom8 = 0.0_rk
   if (oxy > 0) then
     bioom1 = 0.1/secs_pr_day * exp(temp*0.11_rk) * oxy/((0.1_rk*redf(7))+oxy)
     bioom2 = bioom1
     bioom6 = 1.0_rk
   else
     if (no3>0) then
       bioom5 = 5.0_rk
       bioom8 = 1.0_rk
     else
       bioom7 = 1.0_rk
     end if
   end if

!----citical bottom shear stress
        if (tbs.ge.self%BioC(34)) then
!          Rsd=min(self%BioC(35), self%BioC(35) * (tbs**3/(0.1+tbs**3)))
!          Rsd=self%BioC(35)
!          Rsd=min(self%BioC(35), self%BioC(35) * tbs**2 * 500.) ! sets to max=self%BioC(35) when tbs=0.032, else rapid increase from 0. 
          Rsd=min(self%BioC(35), self%BioC(35) * tbs**2 * 100.) ! sets to max=self%BioC(35) when tbs~=0.070, else rapid increase from 0.
          Rds=0.0_rk
        else if (tbs.lt.self%BioC(34)) then
          Rsd=0.0_rk
          Rds=self%BioC(36)
        end if

!---------------------------------------------------------------
!----denitrification parameter in dependence of available oxygen
        if (oxy .gt. 0.0) then
          Rsa=self%BioC(38)*exp(self%BioC(39)*temp)*1.0_rk
          Rsdenit=0.0_rk
        else if (oxy .le. 0.0) then
          Rsdenit=self%BioC(38)*exp(self%BioC(39)*temp)*2.0_rk
          Rsa=0.0_rk
        end if

        !--- sediment 1 total sediment biomass and nitrogen pool
        rhs = Rds*det - Rsd*sed1 - 2.0_rk*Rsa*sed1 - Rsdenit*sed1 &
              -(2.0E-3*self%BioC(37)*sed1)*sed1 !- self%BioC(37)*sed1
        _ADD_BOTTOM_SOURCE_(self%id_sed1, rhs)


        ! oxygen
        flux = -(BioOM6*6.625_rk*2.0_rk*Rsa*sed1 &
                 +BioOM7*6.625_rk*Rsdenit*sed1 &
                 +2.0_rk*BioOM1*Rsa*sed1) &
                *REDF(11)*REDF(16)
        _ADD_BOTTOM_FLUX_(self%id_oxy, flux)

        ! nitrate
        _ADD_BOTTOM_FLUX_(self%id_no3, -BioOM5*Rsdenit*sed1)

        ! detritus
        _ADD_BOTTOM_FLUX_(self%id_det, Rsd*sed1 - Rds*det)

        ! ammonium
        _ADD_BOTTOM_FLUX_(self%id_nh4, (Rsdenit+Rsa)*sed1)

        if (self%couple_co2) then
          _ADD_BOTTOM_FLUX_(self%id_dic, redf(16)*(Rsdenit+2*Rsa)*sed1)
          alk_flux = redf(16)*redf(11)*((Rsdenit+Rsa+bioom5*Rsdenit)*sed1) - 0.5_rk*flux*(1._rk-bioom6)
          _ADD_BOTTOM_FLUX_(self%id_alk, alk_flux)
        end if
        !--try out for phosphate ute 2.6.2010
        Rsa_p=self%BioC(38)*exp(self%BioC(39)*temp)*2.0_rk

        if (oxy.gt.0.0) then
          yt2=oxy/375.0_rk   !normieren des wertes wie in Neumann et al 2002
          yt1=yt2**2.0_rk/(self%BioC(41)**2.0_rk+yt2**2.0_rk)

          _ADD_BOTTOM_FLUX_(self%id_pho,Rsa_p*(1.0_rk-self%BioC(40)*yt1)*sed3)

          !--sed 3 phosphate pool sediment+remineralization-P release
          _ADD_BOTTOM_SOURCE_(self%id_sed3, 2.0_rk*Rsa*sed1 - Rsa_p*(1.0_rk-self%BioC(40)*yt1)*sed3)

        else if (oxy.le.0.0) then
          _ADD_BOTTOM_FLUX_(self%id_pho, Rsa_p*sed3)
          _ADD_BOTTOM_SOURCE_(self%id_sed3, Rsdenit*sed1 - Rsa_p*sed3)
        end if

        ! sediment opal(Si)
!        _ADD_BOTTOM_SOURCE_(self%id_sed2, Rds*opa - Rsd*sed2 - self%BioC(42)*sed2 - ( self%BioC(37)*1000.*(sed2**3/(sed2**3 + 1E+12)) )*sed2)
        _ADD_BOTTOM_SOURCE_(self%id_sed2, 2.0*Rds*opa - Rsd*sed2 - self%BioC(42)*sed2 - (2.0E-3*self%BioC(37)*sed2)*sed2)
        _ADD_BOTTOM_FLUX_(self%id_opa, Rsd*sed2 - 2.0*Rds*opa)
!        _ADD_BOTTOM_FLUX_(self%id_opa, Rsd*sed2 - Rds*opa)
        _ADD_BOTTOM_FLUX_(self%id_sil, self%BioC(42)*sed2)

        if (self%turn_on_additional_diagnostics) then
           _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tbsout, tbs)
        endif
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
   class (type_nersc_ecosmo_operational), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_

   real(rk)                     :: dom,det,diachl,flachl,bgchl
   real(rk)                     :: dia,fla,bg
   real(rk)                     :: my_extinction

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   _GET_(self%id_det, det)
   _GET_(self%id_dom, dom)
   _GET_(self%id_dia, dia)
   _GET_(self%id_fla, fla)
     if (self%use_cyanos) then
   _GET_(self%id_bg, bg)
     else
       bg=0.0_rk
     end if

   my_extinction = self%BioC(4)
   if (self%use_chl) then
     _GET_(self%id_diachl, diachl)
     _GET_(self%id_flachl, flachl)
     if (self%use_cyanos) then
       _GET_(self%id_bgchl, bgchl)
     else
       bgchl=0.0_rk
     end if
     my_extinction = my_extinction + self%BioC(5)*(diachl+flachl+bgchl) + self%extdet*det + self%extdom*dom
   else
     diachl=0.0_rk
     flachl=0.0_rk
     my_extinction = my_extinction + self%BioC(5)*(dia+fla+bg) + self%extdet*det + self%extdom*dom
   end if

   _SET_EXTINCTION_( my_extinction )

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction

! -------------------------------------------------------------------------

   end module fabm_nersc_ecosmo_operational


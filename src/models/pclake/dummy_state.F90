#include "fabm_driver.h"
   module pclake_dummy_state
   
! !DESCRIPTION:
!  this is empty module with all the regiesterd state variables. For the purpose 
!  of exclude some modules but missing certain state variable links 
!
!-------------------------------------------------------------------------------
!BOP
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_dummy_state
!     local state variable identifers
!!    pelagic variables
      type (type_state_variable_id)            :: id_sDIMW,id_sDDetW,id_sNDetW,id_sPDetW,id_sSiDetW
      type (type_state_variable_id)            :: id_sPO4W,id_sPAIMW,id_sNH4W,id_sNO3W,id_sSiO2W,id_sO2W
      type (type_state_variable_id)            :: id_sDDiatW,id_sPDiatW,id_sNDiatW
      type (type_state_variable_id)            :: id_sDGrenW,id_sPGrenW,id_sNGrenW
      type (type_state_variable_id)            :: id_sDBlueW,id_sPBlueW,id_sNBlueW
      type (type_state_variable_id)            :: id_sDZoo,id_sPZoo,id_sNZoo
      type (type_state_variable_id)            :: id_sDFiJv,id_sPFiJv,id_sNFiJv
      type (type_state_variable_id)            :: id_sDFiAd,id_sPFiAd,id_sNFiAd,id_sDPisc 
!!    benthic variables
      type (type_bottom_state_variable_id)     :: id_sPO4S,id_sPAIMS,id_sNH4S,id_sNO3S,id_sDIMS
      type (type_bottom_state_variable_id)     :: id_sDDetS,id_sNDetS,id_sPDetS,id_sSiDetS
      type (type_bottom_state_variable_id)     :: id_sDDiatS,id_sPDiatS,id_sNDiatS
      type (type_bottom_state_variable_id)     :: id_sDGrenS,id_sPGrenS,id_sNGrenS
      type (type_bottom_state_variable_id)     :: id_sDBlueS,id_sPBlueS,id_sNBlueS
      type (type_bottom_state_variable_id)     :: id_sDVeg,id_sNVeg,id_sPVeg
      type (type_bottom_state_variable_id)     :: id_sDBent,id_sPBent,id_sNBent


   contains
!     Model procedures
      procedure :: initialize
   end type type_pclake_dummy_state


!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the inorganic and organic matter model in water
!
! !INTERFACE:

   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
      class (type_pclake_dummy_state), intent(inout), target  :: self
      integer,                          intent(in)            :: configunit
   
!-------------------------------------------------------------------------------
!     LOCAL VARIABLES:
!-------------------------------------------------------------------------------
!!    pelagic
      real(rk)            :: sNH4W_initial
      real(rk)            :: sNO3W_initial
      real(rk)            :: sPO4W_initial
      real(rk)            :: sPAIMW_initial
      real(rk)            :: sSiO2W_initial
      real(rk)            :: sO2W_initial
      real(rk)            :: sDIMW_initial
      real(rk)            :: sDDetW_initial
      real(rk)            :: sPDetW_initial
      real(rk)            :: sNDetW_initial
      real(rk)            :: sSiDetW_initial
      real(rk)            :: sDDiatW_initial
      real(rk)            :: sPDiatW_initial
      real(rk)            :: sNDiatW_initial
      real(rk)            :: sDGrenW_initial
      real(rk)            :: sPGrenW_initial
      real(rk)            :: sNGrenW_initial
      real(rk)            :: sDBlueW_initial
      real(rk)            :: sPBlueW_initial
      real(rk)            :: sNBlueW_initial
      real(rk)            :: sDZoo_initial
      real(rk)            :: sPZoo_initial
      real(rk)            :: sNZoo_initial
      real(rk)            :: sDFiJv_initial
      real(rk)            :: sPFiJv_initial
      real(rk)            :: sNFiJv_initial
      real(rk)            :: sDFiAd_initial
      real(rk)            :: sPFiAd_initial
      real(rk)            :: sNFiAd_initial
      real(rk)            :: sDPisc_initial
!!    benthic
      real(rk)            :: sNH4S_initial
      real(rk)            :: sNO3S_initial
      real(rk)            :: sPO4S_initial
      real(rk)            :: sPAIMS_initial
      real(rk)            :: sDIMS_initial
      real(rk)            :: sDDetS_initial
      real(rk)            :: sPDetS_initial
      real(rk)            :: sNDetS_initial
      real(rk)            :: sSiDetS_initial
      real(rk)            :: sDDiatS_initial
      real(rk)            :: sPDiatS_initial
      real(rk)            :: sNDiatS_initial
      real(rk)            :: sDGrenS_initial
      real(rk)            :: sPGrenS_initial
      real(rk)            :: sNGrenS_initial
      real(rk)            :: sDBlueS_initial
      real(rk)            :: sPBlueS_initial
      real(rk)            :: sNBlueS_initial
      real(rk)            :: sDVeg_initial
      real(rk)            :: sNVeg_initial
      real(rk)            :: sPVeg_initial
      real(rk)            :: sDBent_initial
      real(rk)            :: sPBent_initial
      real(rk)            :: sNBent_initial
      
!-------------------------------------------------------------------------------
!     create namelist
!-------------------------------------------------------------------------------
      namelist /pclake_dummy_state/ sNH4W_initial,sNO3W_initial,sPO4W_initial,sPAIMW_initial,sSiO2W_initial,sO2W_initial, &
                                    & sDIMW_initial,sDDetW_initial,sPDetW_initial,sNDetW_initial,sSiDetW_initial, &
                                    & sDDiatW_initial,sPDiatW_initial,sNDiatW_initial,sDGrenW_initial,sPGrenW_initial,  &
                                    & sNGrenW_initial,sDBlueW_initial,sPBlueW_initial,sNBlueW_initial,&
                                    & sDZoo_initial,sPZoo_initial,sNZoo_initial,sDFiJv_initial,sPFiJv_initial,sNFiJv_initial,& 
                                    & sDFiAd_initial,sPFiAd_initial,sNFiAd_initial,sDPisc_initial, &
!     pelagic above
!     benthic below
                                  & sNH4S_initial,sNO3S_initial,sPO4S_initial,sPAIMS_initial,sDIMS_initial, &
                                  & sDDetS_initial,sNDetS_initial,sPDetS_initial,sSiDetS_initial,&
                                  & sDDiatS_initial,sPDiatS_initial,sNDiatS_initial,sDGrenS_initial,sPGrenS_initial,& 
                                  & sNGrenS_initial,sDBlueS_initial,sPBlueS_initial,sNBlueS_initial, &
                                  & sDVeg_initial,sNVeg_initial,sPVeg_initial,&
                                  & sDBent_initial,sPBent_initial,sNBent_initial
      
      
      
!EOP
!-----------------------------------------------------------------------------
!BOC
!-------------------------------------------------------------------------------
!     initialize the parameters
!-------------------------------------------------------------------------------
!!    pelagic
      sNH4W_initial=0.1_rk
      sNO3W_initial=0.1_rk
      sPO4W_initial=0.01_rk
      sPAIMW_initial=0.0_rk
      sSiO2W_initial=3.0_rk
      sO2W_initial=10.0_rk
      sDIMW_initial=5.0_rk
      sDDetW_initial=2.0_rk
      sPDetW_initial=0.005_rk
      sNDetW_initial=0.05_rk
      sSiDetW_initial=0.02_rk
      sDDiatW_initial=0.5_rk
      sPDiatW_initial=0.005_rk
      sNDiatW_initial=0.05_rk
      sDGrenW_initial=0.5_rk
      sPGrenW_initial=0.005_rk 
      sNGrenW_initial =0.05_rk
      sDBlueW_initial=3.0_rk
      sPBlueW_initial=0.03_rk
      sNBlueW_initial =0.3_rk 
      sDZoo_initial=0.05_rk
      sPZoo_initial=0.0005_rk
      sNZoo_initial=0.0035_rk
      sDFiJv_initial=1.0_rk
      sPFiJv_initial=1.0_rk
      sNFiJv_initial=1.0_rk
      sDFiAd_initial=1.0_rk
      sPFiAd_initial=1.0_rk
      sNFiAd_initial=1.0_rk
      sDPisc_initial=1.0_rk
!!    benthic
      sNH4S_initial=0.02_rk
      sNO3S_initial=0.002_rk
      sPO4S_initial=0.182_rk
      sPAIMS_initial=17.99_rk
      sDIMS_initial=39611.3_rk
      sDDetS_initial=181.7_rk
      sNDetS_initial=4.54_rk
      sPDetS_initial=0.454_rk
      sSiDetS_initial=1.82_rk
      sDDiatS_initial =0.001_rk
      sPDiatS_initial =0.00001_rk
      sNDiatS_initial =0.0001_rk
      sDGrenS_initial =0.001_rk
      sPGrenS_initial =0.00001_rk
      sNGrenS_initial =0.0001_rk
      sDBlueS_initial =0.001_rk
      sPBlueS_initial =0.00001_rk
      sNBlueS_initial =0.0001_rk
      sDVeg_initial=1.0_rk
      sNVeg_initial=0.02_rk
      sPVeg_initial=0.002_rk
      sDBent_initial=1.0_rk
      sPBent_initial=0.1_rk
      sNBent_initial=0.01_rk
      
!-------------------------------------------------------------------------------
!     Read parameters namelist
!-------------------------------------------------------------------------------
      if (configunit>0) read(configunit,nml=pclake_dummy_state,err=99,end=100)
!-------------------------------------------------------------------------------
!     Store parameter values in our own derived type
!-------------------------------------------------------------------------------
!     NB: all rates must be provided in values per day,
!     and are converted here to values per second.

!-----------------------------------------------------------------------------------------------------------------
!     Register local state variable
!-----------------------------------------------------------------------------------------------------------------
!     pelagic variables
      call self%register_state_variable(self%id_sDIMW,'sDIMW','g/m**3','Inorg. matter in water',     &
                                       sDIMW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sDDetW,'sDDetW','g/m**3','detritus dry-weight in water',     &
                                       sDDetW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sPDetW,'sPDetW','g/m**3','detritus phosphrus in water',     &
                                       sPDetW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sNDetW,'sNDetW','g/m**3','detritus nitrogen in water',     &
                                       sNDetW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sSiDetW,'sSiDetW','g/m**3','detritus silica in water',     &
                                       sSiDetW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sPO4W,'sPO4W','g/m**3','Phosphate in water',     &
                                       sPO4W_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sPAIMW,'sPAIMW','g/m**3','Absorbed_P in water',     &
                                       sPAIMW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sNH4W,'sNH4W','g/m**3','Amonia in water',     &
                                       sNH4W_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sNO3W,'sNO3W','g/m**3','Nitrates in water',     &
                                       sNO3W_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sSiO2W,'sSiO2W','g/m**3','Silica dioxide in water',     &
                                       sSiO2W_initial,minimum=_ZERO_,no_river_dilution=.FALSE.)
      call self%register_state_variable(self%id_sO2W,'sO2W','g/m**3','oxygen in water',     &
                                       sO2W_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sDDiatW,'sDDiatW','g/m**3','diatom_D in water',     &
                                       sDDiatW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sPDiatW,'sPDiatW','g/m**3','diatom_P in water',     &
                                       sPDiatW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sNDiatW,'sNDiatW','g/m**3','diatom_N in water',     &
                                       sNDiatW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sDGrenW,'sDGrenW','g/m**3','green_D in water',     &
                                       sDGrenW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sPGrenW,'sPGrenW','g/m**3','green_P in water',     &
                                       sPGrenW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sNGrenW,'sNGrenW','g/m**3','green_N in water',     &
                                       sNGrenW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sDBlueW,'sDBlueW','g/m**3','blue_D in water',     &
                                       sDBlueW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sPBlueW,'sPBlueW','g/m**3','blue_P in water',     &
                                       sPBlueW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sNBlueW,'sNBlueW','g/m**3','blue_N in water',     &
                                    sNBlueW_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sDZoo,'sDZoo','g/m**3','zooplankton DW in water',     &
                                       sDZoo_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sPZoo,'sPZoo','g/m**3','zooplankton P in water',     &
                                       sPZoo_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sNZoo,'sNZoo','g/m**3','zooplankton N in water',     &
                                       sNZoo_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sDFiJv,'sDFiJv','g/m**3','juvenile fish DW in water',     &
                                       sDZoo_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sPFiJv,'sPFiJv','g/m**3','juvenile fish P in water',     &
                                       sPFiJv_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sNFiJv,'sNFiJv','g/m**3','juvenile fish N in water',     &
                                       sNFiJv_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sDFiAd,'sDFiAd','g/m**3','adult fish DW in water',     &
                                       sDFiAd_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sPFiAd,'sPFiAd','g/m**3','adult fish P in water',     &
                                       sPFiAd_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sNFiAd,'sNFiAd','g/m**3','adult fish N in water',     &
                                       sNFiAd_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
      call self%register_state_variable(self%id_sDPisc,'sDPisc','g/m**3','piscivorous fish DW in water',     &
                                       sDPisc_initial,minimum=_ZERO_,no_river_dilution=.TRUE.)
!     benthic
      call self%register_state_variable(self%id_sDIMS,'sDIMS','g/m**2','sediment inorg.Matter',     &
                                       sDIMS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sDDetS,'sDDetS','g/m**2','sediment detritus DW',     &
                                       sDDetS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sNDetS,'sNDetS','g/m**2','sediment detritus N',     &
                                       sNDetS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sPDetS,'sPDetS','g/m**2','sediment detritus P',     &
                                       sPDetS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sSiDetS,'sSiDetS','g/m**2','sediment detritus Si',     &
                                       sSiDetS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sPO4S,'sPO4S','g/m**2','Sediment Phosphate',     &
                                       sPO4S_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sPAIMS,'sPAIMS','g/m**2','SED_Absorbed Phosphate',     &
                                       sPAIMS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sNH4S,'sNH4S','g/m**2','Sediment Amonia',     &
                                       sNH4S_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sNO3S,'sNO3S','g/m**2','Sediment Nitrates',     &
                                       sNO3S_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sDDiatS,'sDDiatS','g/m**2','diatom_D in sediment',     &
                                       sDDiatS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sPDiatS,'sPDiatS','g/m**2','diatom_P in sediment',     &
                                       sPDiatS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sNDiatS,'sNDiatS','g/m**2','diatom_N in sediment',     &
                                       sNDiatS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sDGrenS,'sDGrenS','g/m**2','green_D in sediment',     &
                                       sDGrenS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sPGrenS,'sPGrenS','g/m**2','green_P in sediment',     &
                                       sPGrenS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sNGrenS,'sNGrenS','g/m**2','green_N in sediment',     &
                                       sNGrenS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sDBlueS,'sDBlueS','g/m**2','blue_D in sediment',     &
                                       sDBlueS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sPBlueS,'sPBlueS','g/m**2','blue_P in sediment',     &
                                       sPBlueS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sNBlueS,'sNBlueS','g/m**2','blue_N in sediment',     &
                                       sNBlueS_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sDVeg,'sDVeg','g/m**2','vegetation_dry_weight',     &
                                       sDVeg_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sNVeg,'sNVeg','g/m**2','vegetation_Nitrogen',     &
                                       sNVeg_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sPVeg,'sPVeg','g/m**2','vegetation_Phosphrus',     &
                                       sPVeg_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sDBent,'sDBent','g/m**2','zoobenthos_DW',     &
                                       sDBent_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sPBent,'sPBent','g/m**2','zoobenthos_P',     &
                                       sPBent_initial,minimum=_ZERO_)
      call self%register_state_variable(self%id_sNBent,'sNBent','g/m**2','zoobenthos_N',     &
                                       sNBent_initial,minimum=_ZERO_)
      
      
   return  


99 call self%fatal_error('pclake_dummy_state_init','Error reading namelist pclake_dummy_state')

100 call self%fatal_error('pclake_dummy_state_init','Namelist pclake_dummy_state was not found.')
!

   end subroutine initialize
   
!EOC

!------------------------------------------------------------------------------
   end module pclake_dummy_state

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

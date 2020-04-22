#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_omexdia_p --- Fortran 2003 version of OMEXDIA+P biogeochemical model
!
! !INTERFACE:
   module fabm_hzg_omexdia_p
!
! !DESCRIPTION:
!
! The OMEXDIA+P model is based on the OMEXDIA model (see Soetard et al. 1996a)
! and is intended to simulate early diagenesis in the sea sediments. The major
! difference to the original OMEXDIA is an added phosphorus cycle.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hzg_omexdia_p
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: secs_pr_day = 86400.0_rk
!
! !REVISION HISTORY:!
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hzg_omexdia_p
!     Variable identifiers
      type (type_state_variable_id)        :: id_fdet,id_sdet,id_pdet
      type (type_state_variable_id)        :: id_no3,id_nh3,id_oxy,id_po4,id_odu
      type (type_dependency_id)            :: id_temp
      type (type_diagnostic_variable_id)   :: id_denit,id_adsp

!     Model parameters
      real(rk) :: rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: PCrFdet, PCrSdet, PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox
      real(rk) :: ksO2oduox, ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,kinO2anox

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_hzg_omexdia_p
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the OMEXDIA+P model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the omexdia namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_hzg_omexdia_p),intent(inout),target  :: self
   integer,                   intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
      real(rk) :: rFast, rSlow, NCrFdet, NCrSdet
      real(rk) :: PCrFdet,PCrSdet, PAds, PAdsODU
      real(rk) :: NH3Ads, rnit, ksO2nitri,rODUox
      real(rk) :: ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,kinO2anox
      real(rk) :: fdet_init,sdet_init,oxy_init,odu_init,no3_init,nh3_init
      real(rk) :: pdet_init,po4_init

   namelist /hzg_omexdia_p/  rFast, rSlow, NCrFdet, NCrSdet, &
          PCrFdet, PCrSdet, PAds, PAdsODU, NH3Ads, rnit, ksO2nitri, &
          rODUox,ksO2oduox, ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox, &
          kinO2anox,fdet_init,sdet_init,oxy_init,odu_init,no3_init,nh3_init, &
          pdet_init,po4_init
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Read the namelist
   if (configunit>0) read(configunit,nml=hzg_omexdia_p,err=99,end=100)

   ! Store parameter values in our own derived type
   self%rFast=rFast
   self%rSlow=rSlow
   self%NCrFdet=NCrFdet
   self%NCrSdet=NCrSdet
   self%PCrFdet=PCrFdet
   self%PCrSdet=PCrSdet
   self%PAds=PAds
   self%PAdsODU=PAdsODU
   self%NH3Ads=NH3Ads
   self%rnit=rnit
   self%ksO2nitri=ksO2nitri
   self%rODUox=rODUox
   self%ksO2oduox=ksO2oduox
   self%ksO2oxic=ksO2oxic
   self%ksNO3denit=ksNO3denit
   self%kinO2denit=kinO2denit
   self%kinNO3anox=kinNO3anox
   self%kinO2anox=kinO2anox

   ! Register state variables
   call self%register_state_variable(self%id_fdet,'fdet','mmolC/m**3','fast detritus C',     &
                                    fdet_init,minimum=0.0_rk)
   call self%set_variable_property(self%id_fdet,'particulate',.true.)

   call self%register_state_variable(self%id_sdet,'sdet','mmolC/m**3','slow detritus C', &
                                    sdet_init,minimum=0.0_rk)
   call self%set_variable_property(self%id_sdet,'particulate',.true.)

   call self%register_state_variable(self%id_pdet,'pdet','mmolP/m**3','detritus-P',     &
                                    pdet_init,minimum=0.0_rk)
   call self%set_variable_property(self%id_pdet,'particulate',.true.)

   call self%register_state_variable(self%id_po4,'po4','mmolP/m**3','dissolved phosphate', &
                                    po4_init,minimum=0.0_rk, &
                                    standard_variable=standard_variables%mole_concentration_of_phosphate)
   call self%set_variable_property(self%id_po4,'particulate',.false.)

   call self%register_state_variable(self%id_no3,'no3','mmolN/m**3','dissolved nitrate',     &
                                    no3_init,minimum=0.0_rk, &
                                    standard_variable=standard_variables%mole_concentration_of_nitrate)
   call self%set_variable_property(self%id_no3,'particulate',.false.)

   call self%register_state_variable(self%id_nh3,'nh3','mmolN/m**3','dissolved ammonium', &
                                    nh3_init,minimum=0.0_rk, &
                                    standard_variable=standard_variables%mole_concentration_of_ammonium)
   call self%set_variable_property(self%id_nh3,'particulate',.false.)

   call self%register_state_variable(self%id_oxy,'oxy','mmolO2/m**3','dissolved oxygen',     &
                                    oxy_init,minimum=0.0_rk)
   call self%set_variable_property(self%id_oxy,'particulate',.false.)

   call self%register_state_variable(self%id_odu,'odu','mmol/m**3','dissolved reduced substances', &
                                    odu_init,minimum=0.0_rk)
   call self%set_variable_property(self%id_odu,'particulate',.false.)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_adsp,'adsP','mmolP/m**3', &
         'phosphate adsorption', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_denit,'denit','mmol/m**3/d', &
         'denitrification rate', output=output_instantaneous)

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

   return

99 call self%fatal_error('hzg_omexdia_p_initialize','Error reading namelist hzg_omexdia_p.')

100 call self%fatal_error('hzg_omexdia_p_initialize','Namelist hzg_omexdia_p was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of OMEXDIA+P model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_hzg_omexdia_p),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister & Kai Wirtz
!
! !LOCAL VARIABLES:
   real(rk) :: fdet,sdet,oxy,odu,no3,nh3,pdet,po4
   real(rk) :: temp_celsius,temp_kelvin,f_T,E_a
   real(rk) :: radsP,Oxicminlim,Denitrilim,Anoxiclim,Rescale,rP
   real(rk),parameter :: relaxO2=0.04_rk
   real(rk),parameter :: T0 = 288.15_rk ! reference Temperature fixed to 15 degC
   real(rk),parameter :: Q10b = 1.5_rk
   real(rk) :: CprodF,CprodS,Cprod,Nprod,Pprod
   real(rk) :: AnoxicMin,Denitrific,OxicMin,Nitri,OduDepo,OduOx,pDepo
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_temp,temp_celsius)
   _GET_(self%id_fdet,fdet)
   _GET_(self%id_sdet,sdet)
   _GET_(self%id_pdet,pdet)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_odu,odu)
   _GET_(self%id_no3,no3)
   _GET_(self%id_nh3,nh3)
   _GET_(self%id_po4,po4)

   temp_kelvin = 273.15_rk + temp_celsius
   E_a=0.1_rk*log(Q10b)*T0*(T0+10.0_rk);
   f_T = 1.0_rk*exp(-E_a*(1.0_rk/temp_kelvin - 1.0_rk/T0))

   Oxicminlim = oxy/(oxy+self%ksO2oxic+relaxO2*(nh3+odu))                ! limitation terms
   Denitrilim = (1.0_rk-oxy/(oxy+self%kinO2denit)) * NO3/(no3+self%ksNO3denit)
   Anoxiclim  = (1.0_rk-oxy/(oxy+self%kinO2anox)) * (1.0_rk-no3/(no3+self%kinNO3anox))
   Rescale    = 1.0_rk/(Oxicminlim+Denitrilim+Anoxiclim)

   CprodF = self%rFast * fdet
   CprodS = self%rSlow * sdet
   Cprod  = CprodF + CprodS
   Nprod  = CprodF * self%NCrFdet + CprodS * self%NCrSdet


! PO4-adsorption ceases when critical capacity is reached
! [FeS] approximated by ODU
   radsP  = self%PAds * self%rSlow * (po4*max(odu,self%PAdsODU))
   rP    = self%rFast * (1.0_rk - Oxicminlim)
   Pprod  = rP * pdet

! Oxic mineralisation, denitrification, anoxic mineralisation
! then the mineralisation rates
   OxicMin    = Cprod*Oxicminlim*Rescale        ! oxic mineralisation
   Denitrific = Cprod*Denitrilim*Rescale        ! Denitrification
   AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

! reoxidation and ODU deposition
   Nitri      = f_T * self%rnit   * nh3 * oxy/(oxy + self%ksO2nitri + relaxO2*(fdet + odu))
   OduOx      = f_T * self%rODUox * odu * oxy/(oxy + self%ksO2oduox + relaxO2*(nh3 + fdet))

!  pDepo      = min(1.0_rk,0.233_rk*(wDepo)**0.336_rk )
   pDepo      = 0.0_rk
   OduDepo    = AnoxicMin*pDepo

#define _CONV_UNIT_ /secs_pr_day
! reaction rates
   _ADD_SOURCE_(self%id_fdet, -f_T * CprodF _CONV_UNIT_)
   _ADD_SOURCE_(self%id_sdet, -f_T * CprodS _CONV_UNIT_)
   _ADD_SOURCE_(self%id_oxy , (-OxicMin - 2.0_rk* Nitri - OduOx) _CONV_UNIT_) !RH 1.0->150/106*OxicMin (if [oxy]=mmolO2/m**3)
   _ADD_SOURCE_(self%id_no3 , (-0.8_rk*Denitrific + Nitri) _CONV_UNIT_)     !RH 0.8-> ~104/106?
   _ADD_SOURCE_(self%id_nh3 , (f_T * Nprod - Nitri) / (1.0_rk + self%NH3Ads) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_odu , (AnoxicMin - OduOx - OduDepo) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_po4 , (f_T * Pprod - radsP) _CONV_UNIT_)
   _ADD_SOURCE_(self%id_pdet, (radsP - f_T * Pprod) _CONV_UNIT_)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_denit,Denitrific)
   _SET_DIAGNOSTIC_(self%id_adsp ,radsP)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

   end module fabm_hzg_omexdia_p


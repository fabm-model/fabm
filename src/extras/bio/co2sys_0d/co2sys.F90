!$Id: bio_npzd_0d.F90,v 1.6 2009-05-10 18:36:38 jorn Exp $
#include "rmbm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_co2_sys_0d --- 0D CO2/carbonate system model
!
! !INTERFACE:
module rmbm_co2sys
!
! !DESCRIPTION:
! CO2 system model based on PML code.
!
! !USES:
   use rmbm_types
   use rmbm_driver

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_co2_sys, init_bio_co2_sys_0d, do_bio_co2_sys_0d, update_air_sea_co2_sys_0d
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!
! !PUBLIC DERIVED TYPES:
   type type_co2_sys
      integer :: id_dic
      integer :: id_temp, id_salt, id_pres, id_wnd
      integer :: id_ph, id_TA, id_pco2, id_cco2, id_CarbA, id_Bicarb, id_Carb, id_Om_cal, id_Om_arg
      REALTYPE  :: TA_offset, TA_slope, pCO2a
   end type
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio_co2_sys_0d(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the bio namelist is read and 
!  various variables (rates and settling velocities) 
!  are transformed into SI units.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_co2_sys), intent(out)     :: self
   type (type_model_info),intent(inout) :: modelinfo
   integer,               intent(in )   :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE :: dic_initial
   REALTYPE  :: TA_offset = 520.1, TA_slope = 51.24, pCO2a = 390.
   namelist /bio_co2_sys_nml/ dic_initial, TA_offset, TA_slope, pCO2a
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=bio_co2_sys_nml,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day, and are converted here to values per second.
   self%TA_offset = TA_offset
   self%TA_slope  = TA_slope
   self%pCO2a     = pCO2a
      
   self%id_dic = register_state_variable(modelinfo,'DIC','mol C/m**3','total dissolved inorganic carbon', &
                                    dic_initial,minimum=_ZERO_)
                                    
   self%id_ph     = register_diagnostic_variable(modelinfo, 'pH',    '-',        'pH',                           2)
   self%id_TA     = register_diagnostic_variable(modelinfo, 'TA',    'mmol/m**3','alkalinity',                   2)
   self%id_pco2   = register_diagnostic_variable(modelinfo, 'pCO2',  'ppm',      'CO2 partial pressure',         2)
   self%id_cco2   = register_diagnostic_variable(modelinfo, 'cCO2',  'mmol/m**3','CO2 concentration',            2)
   self%id_CarbA  = register_diagnostic_variable(modelinfo, 'CarbA', 'mmol/m**3','carbonic acid concentration',  2)
   self%id_Bicarb = register_diagnostic_variable(modelinfo, 'Bicarb','mmol/m**3','bicarbonate ion concentration',2)
   self%id_Carb   = register_diagnostic_variable(modelinfo, 'Carb',  'mmol/m**3','carbonate ion concentration',  2)
   self%id_Om_cal = register_diagnostic_variable(modelinfo, 'Om_cal','-',        'calcite saturation state',     2)
   self%id_Om_arg = register_diagnostic_variable(modelinfo, 'Om_arg','-',        'aragonite saturation state',   2)

   return

99 call fatal_error('init_bio_co2_sys_0d','I could not read namelist bio_co2_sys_nml')
   
   end subroutine init_bio_co2_sys_0d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of carbonate system model
!
! !INTERFACE:
   subroutine do_bio_co2_sys_0d(self,state,environment,LOCATION,dy,diag)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_co2_sys),    intent(in) :: self
   type (type_state),      intent(in) :: state(:)
   type (type_environment),intent(in) :: environment
   LOCATIONTYPE                       :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout),dimension(:) :: dy,diag
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   ! Environment
   REALTYPE :: temp, salt, pres
   REALTYPE :: dic, TA
   REALTYPE :: PCO2WATER, pH, HENRY, PCO2X, ca, bc, cb,Om_cal,Om_arg
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Get current value for total dissolved inorganic carbon (our own state variable).
   dic = state(self%id_dic)%data(LOCATION)

   ! Get environmental variables.
   temp = environment%temp(LOCATION)
   salt = environment%salt(LOCATION)
   pres = environment%pres(LOCATION)

   ! Linearly approximate alkalinity from salinity.
   TA = self%TA_offset + self%TA_slope*salt

   ! Calculate carbonate system equilibrium.
   call CO2DYN(dic/1.0D6, TA/1.0D6, temp, salt, PCO2WATER, pH, HENRY, PCO2X, ca, bc, cb)

   !Call carbonate saturation state subroutine to calculate calcite and aragonite calcification states.
   call CaCO3_Saturation (temp, salt, pres*10.d0, cb, Om_cal, Om_arg)
   
   ! Store diagnostic variables.
   diag(self%id_ph)     = ph
   diag(self%id_TA)     = TA
   diag(self%id_pco2)   = PCO2WATER*1.0D6   ! To ppm
   diag(self%id_cco2)   = PCO2X             ! Already in mmol/m**3
   diag(self%id_CarbA)  = ca       *1.0D6   ! To mmol/m**3
   diag(self%id_Bicarb) = bc       *1.0D6   ! To mmol/m**3
   diag(self%id_Carb)   = cb       *1.0D6   ! To mmol/m**3
   diag(self%id_Om_cal) = Om_cal
   diag(self%id_Om_arg) = Om_arg

   end subroutine do_bio_co2_sys_0d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Air-sea exchange for the carbonate system model
!
! !INTERFACE:
   subroutine update_air_sea_co2_sys_0d(self,state,environment,LOCATION,numc,flux)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_co2_sys), intent(in)      :: self
   type (type_state),      intent(in),pointer :: state(:)
   type (type_environment),intent(in),pointer :: environment
   LOCATIONTYPE                         :: LOCATION
   integer,intent(in)                   :: numc
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, intent(inout)              :: flux(1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   ! Environment
   REALTYPE :: temp, salt, wnd
   
   ! State
   REALTYPE :: dic, TA
   
   ! Temporary variables
   REALTYPE :: PCO2WATER, pH, HENRY, PCO2X, ca, bc, cb, fl

   ! Parameters
   REALTYPE, parameter :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   temp = environment%temp(LOCATION)
   salt = environment%salt(LOCATION)
   wnd  = environment%wind_sf(LOCATION2D)

   dic = state(self%id_dic)%data(LOCATION)

   ! Linearly approximate alkalinity from salinity.
   TA = self%TA_offset + self%TA_slope*salt

   ! Calculate carbonate system equilibrium.
   call CO2DYN(dic/1.0D6, TA/1.0D6, temp, salt, PCO2WATER, pH, HENRY, PCO2X, ca, bc, cb)

   call Air_sea_exchange(temp, wnd, PCO2WATER*1.0D6, self%pCO2a, Henry, fl)
   
   flux(self%id_dic) = fl/secs_pr_day

   end subroutine update_air_sea_co2_sys_0d
!EOC

!-----------------------------------------------------------------------

end module rmbm_co2sys

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

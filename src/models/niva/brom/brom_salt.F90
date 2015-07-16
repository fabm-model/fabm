#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_salt
!
! !DESCRIPTION:
!
! !USES:

   use fabm_types

   implicit none

!  default: all is private.
   private

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_salt
!     Variable identifiers
      type (type_state_variable_id)        :: id_Na,id_Cl,id_NaCl
      type (type_dependency_id)            :: id_temp
            
 real(rk) :: K_NaCl != 37.19e12 ! NaCl conditional equilibrium constant %  37.19 mol2/l2    (Krumgalz, Millero, 1989 )   !!1.8e-11 (M) (Internet)            
 real(rk) :: K_NaCl_diss !=2. ! Specific rate of dissolution of NaCl kg/m3/d !!!(1/day)   !2.5 X 10-1 yr-1 (vanCap06)
 real(rk) :: K_NaCl_form !=1. ! Specific rate of formation of NaCl kg/m3/d !!!(1/day)  !! 1. X 10-4 yr-1(vanCap06)    
 real(rk) :: Wsed != 5. !1Rate of sinking of detritus (POP, PON)d-1 !!  Wdetr=1.5 (Savchuk, Wulff,1996),!Wdetr= 3.5; 20. (Gregoire,2000)

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
   class (type_niva_brom_salt), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%register_state_variable(self%id_Na, 'Na', 'mmol/m**3','sodium', minimum=0.0_rk)
   call self%register_state_variable(self%id_Cl, 'Cl', 'mmol/m**3','chloride', minimum=0.0_rk)
   call self%register_state_variable(self%id_NaCl, 'NaCl', 'mmol/m**3','sodium chloride', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)

   call self%get_parameter(self%K_NaCl,'K_NaCl','umol2/l2',' NaCl conditional equilibrium constant',default=37.19e12_rk)
   call self%get_parameter(self%K_NaCl_diss,'K_NaCl_diss','umol/l/d',' Specific rate of dissolution of NaCl',default=1._rk)
   call self%get_parameter(self%K_NaCl_form,'K_NaCl_form','umol/l/d',' Specific rate of formation of NaCl',default=1._rk)
   call self%get_parameter(self%K_NaCl_form,'Wsed','m/d',' Specific rate of formation of NaCl',default=1._rk)
   
 !  real(rk) :: K_NaCl != 37.19e12 ! NaCl conditional equilibrium constant %  37.19 mol2/l2    (Krumgalz, Millero, 1989 )   !!1.8e-11 (M) (Internet)            
 !real(rk) :: K_NaCl_diss !=2. ! Specific rate of dissolution of NaCl kg/m3/d !!!(1/day)   !2.5 X 10-1 yr-1 (vanCap06)
 !real(rk) :: K_NaCl_form !=1. ! Specific rate of formation of NaCl kg/m3/d !!!(1/day)  !! 1. X 10-4 yr-1(vanCap06)    
 !real(rk) :: Wsed != 5. !1Rate of sinking of detritus (POP, PON)d-1 !!  Wdetr=1.5 (Savchuk, Wulff,1996),!Wdetr= 3.5; 20. (Gregoire,2000)

   
   
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   
    self%dt = 86400.

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
   class (type_niva_brom_salt),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:
   real(rk) :: temp,Na,Cl,NaCl
   real(rk) :: Om_NaCl,NaCl_prec,NaCl_diss
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
   _GET_(self%id_temp,temp)              ! temperature

   _GET_(self%id_Na,Na)
   _GET_(self%id_Cl,Cl)
   _GET_(self%id_NaCl,NaCl)
   
   
! NaCl precipitation/dissolution
        Om_NaCl=Na*Cl/(self%K_NaCl)
!% Ni2+ + 2HCO3- = MnCO3 +CO2 +H2O (vanCap)
        NaCl_prec=self%K_NaCl_form*max(0._rk,(Om_NaCl-1._rk))
!% 
        NaCl_diss=self%K_NaCl_diss*max(0._rk,(1._rk-Om_NaCl))   
   
   
   _SET_ODE_(self%id_Na,-NaCl_prec+NaCl_diss)
   _SET_ODE_(self%id_Cl,-NaCl_prec+NaCl_diss)
   _SET_ODE_(self%id_NaCl,+NaCl_prec-NaCl_diss)
  
! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

end module
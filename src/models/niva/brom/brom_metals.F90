#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_metals
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
!  Original author(s): Evgeniy Yakushev, Svetlana Pakhomova, Jorn Bruggeman
!

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_metals
!     Variable identifiers
      type (type_state_variable_id)        :: id_Cu,id_CuS,id_Cu_prt
      type (type_dependency_id)            :: id_temp
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
   class (type_niva_brom_metals), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%register_state_variable(self%id_Cu, 'Cu', 'mol/m**3','copper', minimum=0.0_rk)
   call self%register_state_variable(self%id_CuS, 'Cl', 'mol/m**3','copper sulphide', minimum=0.0_rk)
   call self%register_state_variable(self%id_Cu_prt, 'Cu_prt', 'mol/m**3','copper in complexes with POM', minimum=0.0_rk)

   call self%register_dependency(self%id_temp,standard_variables%temperature)

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
   class (type_niva_brom_metals),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:
   real(rk) :: temp,Cu, CuS, Cu_prt
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
   _GET_(self%id_temp,temp)              ! temperature

   !_GET_(self%id_Cu,Cu) 
   !_GET_(self%id_CuS,Cl) 
   !_GET_(self%id_Cu_prt,Cu_prt)
   
   
Cu=Cu+CuS*0. !+H2S
   
   _SET_ODE_(self%id_Cu,0.0_rk)
   _SET_ODE_(self%id_CuS,0.0_rk)
   _SET_ODE_(self%id_Cu_prt,0.0_rk)
  
! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

end module
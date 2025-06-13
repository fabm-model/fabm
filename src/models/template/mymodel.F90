#include "fabm_driver.h"

module template_mymodel

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_template_mymodel
      ! Add variable identifiers and parameters here.
   contains
      procedure :: initialize
      ! Reference model procedures here.
   end type

contains

   subroutine initialize(self, configunit)
      class (type_template_mymodel), intent(inout), target :: self
      integer,                       intent(in)            :: configunit
 
      ! Register model parameters and variables here.
   end subroutine initialize

   ! Add model subroutines here.

end module

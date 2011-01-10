!$Id: rmbm_driver.F90 91 2010-12-10 11:59:58Z jornbr $
#include "cppdefs.h"

module rmbm_driver

   implicit none

   public fatal_error,log_message

   private
   
contains

   subroutine fatal_error(routine,errormsg)
     character(len=*), intent(in) :: routine,errormsg
     
     FATAL trim(routine)//': '//trim(errormsg)
     stop 1
   end subroutine fatal_error

   subroutine log_message(msg)
     character(len=*), intent(in) :: msg
     
     write (*,*) trim(msg)
   end subroutine log_message

end module rmbm_driver
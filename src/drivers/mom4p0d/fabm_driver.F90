#include "fabm_driver.h"

module fabm_driver

   use mpp_mod,    only: mpp_error, FATAL, stdout

   implicit none

   public fatal_error,log_message

   private

contains

   subroutine fatal_error(routine,errormsg)
     character(len=*), intent(in) :: routine,errormsg

     call mpp_error(FATAL, trim(routine)//': '//trim(errormsg))
   end subroutine fatal_error

   subroutine log_message(msg)
     character(len=*), intent(in) :: msg

     write (stdout(),*) trim(msg)
   end subroutine log_message

end module fabm_driver
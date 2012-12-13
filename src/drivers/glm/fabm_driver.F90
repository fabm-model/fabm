!$Id: fabm_driver.F90 91 2010-12-10 11:59:58Z jornbr $
#include "cppdefs.h"

MODULE fabm_driver

   IMPLICIT NONE

   PUBLIC fatal_error,log_message

   PRIVATE

CONTAINS

SUBROUTINE fatal_error(routine,errormsg)
   CHARACTER(len=*),INTENT(in) :: routine,errormsg

   FATAL TRIM(routine)//': '//TRIM(errormsg)
   STOP 1
END SUBROUTINE fatal_error

SUBROUTINE log_message(msg)
   CHARACTER(len=*),INTENT(in) :: msg

   write (*,*) TRIM(msg)
END SUBROUTINE log_message

END MODULE fabm_driver

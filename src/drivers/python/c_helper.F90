module fabm_c_helper

   use iso_c_binding, only: c_char, C_NULL_CHAR

   implicit none

contains

   subroutine copy_to_c_string(string,cstring)
      character(len=*),      intent(in)  :: string
      character(kind=c_char),intent(out) :: cstring(:)
      integer i,n
      n = min(len_trim(string),size(cstring)-1)
      do i=1,n
         cstring(i) = string(i:i)
      end do
      cstring(n+1) = C_NULL_CHAR
   end subroutine

   function logical2int(value) result(ivalue)
      logical,intent(in) :: value
      integer            :: ivalue
      if (value) then
         ivalue = 1
      else
         ivalue = 0
      end if
   end function

   function int2logical(ivalue) result(value)
      integer,intent(in) :: ivalue
      logical            :: value
      value = ivalue/=0
   end function

end module
module shared
   use fabm
   use fabm_types

   implicit none

   public

   type (type_model),pointer :: model

   real(rk),allocatable :: cc(:)

   real(rk),target :: temp,salt,par
   real(rk)        :: latitude,longitude
   character(len=80) :: title

end module
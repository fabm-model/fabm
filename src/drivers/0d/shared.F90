module shared
   use fabm
   use fabm_types

   implicit none

   public

   type (type_model),pointer :: model

   real(rk),allocatable :: cc(:)

   real(rk),target :: temp,salt,par
   real(rk)        :: latitude,longitude,column_depth
   character(len=80) :: title

   type type_input_data
      character(len=attribute_length) :: variable_name = ''
      real(rk)                        :: value         = 0.0_rk
      integer                         :: ncid
      type (type_input_data),pointer  :: next          => null()
   end type
   type (type_input_data), pointer, save :: first_input_data => null()

end module
module shared
   use fabm,only: type_model
   use fabm_types,only: attribute_length,rk

   implicit none

   public

   type (type_model),pointer :: model

   real(rk),allocatable,target :: cc(:)

   real(rk),target :: temp,salt,par
   real(rk)        :: latitude,longitude,column_depth
   character(len=80) :: title

   type type_input_data
      character(len=attribute_length) :: variable_name = ''
      real(rk)                        :: value         = 0.0_rk
      integer                         :: ncid          = -1
      type (type_input_data),pointer  :: next          => null()
   end type
   type (type_input_data), pointer, save :: first_input_data => null()

   logical :: compute_conserved_quantities
   real(rk), allocatable, dimension(:) :: totals,int_change_in_totals

end module

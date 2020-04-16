module shared
   use fabm, only: type_fabm_model
   use fabm_types, only: attribute_length, rk
   use input, only: type_scalar_input, type_scalar_input_list

   implicit none

   public

   class (type_fabm_model), save, pointer :: model => null()

   real(rk), allocatable, target :: cc(:)

   type (type_scalar_input), save, target :: light, temp, salt

   type (type_scalar_input_list), save :: extra_inputs

   real(rk), save, target :: latitude, longitude, column_depth, par
   character(len=80) :: title

   logical :: compute_conserved_quantities
   real(rk), allocatable, dimension(:) :: totals,int_change_in_totals

end module

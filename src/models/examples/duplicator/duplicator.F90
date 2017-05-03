#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: examples_duplicator
!
! !INTERFACE:
   module examples_duplicator
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_examples_duplicator
   contains
      procedure :: initialize
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the nutrient componet
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the examples_npzd_nut namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_examples_duplicator), intent(inout), target :: self
   integer,                          intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   character(len=256) :: model,configfile,parameter
   integer            :: number
   real(rk)           :: minimum,maximum,value
   namelist /examples_duplicator/ model,number,configfile,parameter,minimum,maximum

   integer :: i
   character(len=64) :: instancename
   class (type_base_model),pointer :: childmodel
   integer :: childunit
!EOP
!-----------------------------------------------------------------------
!BOC
   model = ''
   configfile = ''
   number = 0
   parameter = ''
   minimum = 0.0_rk
   maximum = 0.0_rk

   ! Read the namelist
   if (configunit>0) read(configunit,nml=examples_duplicator,err=99)

   if (model=='')        call self%fatal_error('initialize','parameter "model" must be provided.')
   if (number==0)        call self%fatal_error('initialize','parameter "number" must be provided.')
   if (configfile=='')   call self%fatal_error('initialize','parameter "configfile" must be provided.')
   if (maximum<=minimum) call self%fatal_error('initialize','parameters "minimum" and "maximum" must be provided, &
                                                            &and "maximum" must exceed "minimum".')

   childunit = get_free_unit()

   do i=1,number
      ! Create a unique instance name for the child model.
      write (unit=instancename,fmt='(a,i2.2)') trim(model),i

      ! Create the child model.
      call factory%create(model,childmodel)
      if (.not.associated(childmodel)) call self%fatal_error('initialize','Model class named "'//model//'" was not found.')

      ! Generate a parameter value and transfer it to the child model.
      call random_number(value)
      value = minimum + value*(maximum-minimum)
      call childmodel%parameters%set_real(parameter,value)

      ! Initialize the child model.
      open(unit=childunit,file=configfile,action='read',status='old')
      call self%add_child(childmodel,instancename,configunit=childunit)
      close(childunit)

      if (.not.childmodel%parameters%retrieved%contains(parameter)) &
         call self%fatal_error('initialize','Parameter "'//trim(parameter)//'" was not registered by model "'//trim(model)//'".')
   end do

   return

99 call self%fatal_error('initialize','Error reading namelist examples_duplicator')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------

   end module examples_duplicator

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#ifdef _FABM_F2003_

#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_examples_duplicator
!
! !INTERFACE:
   module fabm_examples_duplicator
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types
   use fabm_driver
   
   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_examples_duplicator
   
      character(len=256) :: parameter
      real(rk)           :: minimum,maximum
      logical            :: found
      
      contains
      
      procedure :: initialize
      procedure :: get_real_parameter
      
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
   real(rk)           :: minimum,maximum
   namelist /examples_duplicator/ model,number,configfile,parameter,minimum,maximum
   
   integer :: i
   character(len=64) :: instancename
   integer, parameter :: childunit = 10000
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
   read(configunit,nml=examples_duplicator,err=99)
   
   if (model=='') call fatal_error('fabm_examples_duplicator::initialize','parameter "model" must be provided.')
   if (number==0) call fatal_error('fabm_examples_duplicator::initialize','parameter "number" must be provided.')
   if (configfile=='') call fatal_error('fabm_examples_duplicator::initialize','parameter "configfile" must be provided.')
   if (maximum<=minimum) call fatal_error('fabm_examples_duplicator::initialize','parameters "minimum" and "maximum" must be provided, and "maximum" must exceed "minimum".')
   
   self%parameter = parameter
   self%minimum = minimum
   self%maximum = maximum
   self%found = .false.
   
   do i=1,number
      ! Create a unique instance name for the child model.
      write (unit=instancename,fmt='(a,i2.2)') trim(model),i
      
      ! Create and initialize the child model.
      open(unit=childunit,file=configfile,action='read',status='old')
      call self%add_child(model,childunit,instancename)
      close(childunit)
   end do
   
   if (.not.self%found) call fatal_error('fabm_examples_duplicator::initialize','Parameter "'//trim(self%parameter)//'" was not registered by model "'//trim(model)//'".')

   return

99 call fatal_error('fabm_examples_duplicator::initialize','Error reading namelist examples_duplicator')

   end subroutine initialize
!EOC

   recursive subroutine get_real_parameter(self,value,name,units,long_name,scale_factor,default,path)
   ! !INPUT PARAMETERS:
      class (type_examples_duplicator), intent(inout), target :: self
      real(rk),        intent(inout)       :: value
      character(len=*),intent(in)          :: name
      character(len=*),intent(in),optional :: units,long_name,path
      real(rk),        intent(in),optional :: scale_factor,default
      
      call self%type_base_model%get_real_parameter(value,name,units,long_name,scale_factor,default,path)
      
      if (name==self%parameter) then
         call random_number(value)
         value = self%minimum + value*(self%maximum-self%minimum)
         self%found = .true.
      end if
   end subroutine

!-----------------------------------------------------------------------

   end module fabm_examples_duplicator

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#endif
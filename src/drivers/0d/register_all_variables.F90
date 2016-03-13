#ifdef NETCDF4
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: register_all_variables
!
! !INTERFACE:
   module register_all_variables
!
! !DESCRIPTION:
!
! !USES:
   use field_manager
   use fabm
   use fabm_types, only: rk,output_none
   use fabm_properties

   IMPLICIT NONE
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: do_register_all_variables
!
! !PUBLIC DATA MEMBERS:
   type (type_field_manager), public, target :: fm
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !PRIVATE DATA MEMBERS
   integer :: N=0
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_register_all_variables
!
! !INTERFACE:
   subroutine do_register_all_variables(lat,lon,par,temp,salt,model)
!
! !USES:
   IMPLICIT NONE
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: lat,lon
   real(rk)                            :: par,temp,salt
   class (type_model), intent(inout)   :: model
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Jorn Bruggeman
!
! !LOCAL VARIABLES:
!EOP
!-------------------------------------------------------------------------
!BOC
   LEVEL2 'register_all_variables'
   call register_coordinate_variables(lat,lon)
   call register_environment_variables(par,temp,salt)
   call register_fabm_variables(model)
!KB   call fm%list
!   LEVEL2 'registrated ',N,'variables'
   return
   end subroutine do_register_all_variables
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Coordinate variable registration 
!
! !INTERFACE:
   subroutine register_coordinate_variables(lat,lon)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(rk), intent(in) :: lat,lon
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL3 'register_coordinate_variables'

!  register - dimension
   call fm%register_dimension('lon',1,id=id_dim_lon)
   call fm%register_dimension('lat',1,id=id_dim_lat)
   call fm%register_dimension('time',id=id_dim_time)
   call fm%initialize(prepend_by_default=(/id_dim_lon,id_dim_lat/),append_by_default=(/id_dim_time/))
   call fm%register('lon','degrees_east','longitude',dimensions=(/id_dim_lon/),no_default_dimensions=.true.,data0d=lon,coordinate_dimension=id_dim_lon)
   call fm%register('lat','degrees_north','latitude',dimensions=(/id_dim_lat/),no_default_dimensions=.true.,data0d=lat,coordinate_dimension=id_dim_lat)

   return
   end subroutine register_coordinate_variables
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: airsea variable registration 
!
! !INTERFACE:
   subroutine register_environment_variables(par,temp,salt)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(rk)        :: par,temp,salt
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL3 'register_environment_variables'
   call fm%register('par','W/m^2','par',standard_name='downwelling_photosynthetic_radiative_flux',data0d=par)
   call fm%register('temp','Celsius','temperature',standard_name='sea_water_temperature',data0d=temp)
   call fm%register('salt','1e-3','salinity',standard_name='sea_water_practical_salinity',data0d=salt)
   return
   end subroutine register_environment_variables
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: airsea variable registration 
!
! !INTERFACE:
   subroutine register_fabm_variables(model)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   class (type_model), intent(inout)   :: model
!
! !LOCAL VARIABLES:
   integer          :: i
   logical          :: in_output
   real(rk),pointer :: pdata
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL3 'register_fabm_variables'

   ! state variables
   do i=1,size(model%state_variables)
      in_output = register(model%state_variables(i))
   end do
   do i=1,size(model%bottom_state_variables)
      in_output = register(model%bottom_state_variables(i))
   end do
   do i=1,size(model%surface_state_variables)
      in_output = register(model%surface_state_variables(i))
   end do

   ! diagnostic variables
   do i=1,size(model%diagnostic_variables)
      in_output = register(model%diagnostic_variables(i))
      if (in_output) model%diagnostic_variables(i)%save = .true.
   end do
   do i=1,size(model%horizontal_diagnostic_variables)
      in_output = register(model%horizontal_diagnostic_variables(i))
      if (in_output) model%horizontal_diagnostic_variables(i)%save = .true.
   end do

   do i=1,size(model%state_variables)
      pdata => model%get_data(model%state_variables(i)%globalid)
      call fm%send_data(model%state_variables(i)%name, pdata)
   end do
   do i=1,size(model%bottom_state_variables)
      pdata => model%get_data(model%bottom_state_variables(i)%globalid)
      call fm%send_data(model%bottom_state_variables(i)%name, pdata)
   end do
   do i=1,size(model%surface_state_variables)
      pdata => model%get_data(model%surface_state_variables(i)%globalid)
      call fm%send_data(model%surface_state_variables(i)%name, pdata)
   end do

   do i=1,size(model%diagnostic_variables)
      if (model%diagnostic_variables(i)%save) &
         call fm%send_data(model%diagnostic_variables(i)%name, fabm_get_bulk_diagnostic_data(model,i))
   end do
   do i=1,size(model%horizontal_diagnostic_variables)
      if (model%horizontal_diagnostic_variables(i)%save) &
         call fm%send_data(model%horizontal_diagnostic_variables(i)%name, fabm_get_horizontal_diagnostic_data(model,i))
   end do

   contains

   function register(variable) result(used)
      class (type_external_variable),intent(in) :: variable
      logical                                   :: used

      integer                        :: output_level
      type (type_field),     pointer :: field
      class (type_property), pointer :: property

      output_level = output_level_default
      if (variable%output==output_none) output_level = output_level_debug
      call fm%register(variable%name, variable%units, variable%long_name, &
                       minimum=variable%minimum, maximum=variable%maximum, fill_value=variable%missing_value, &
                       category='fabm'//variable%target%owner%get_path(), output_level=output_level, used=used, field=field)
      property => variable%properties%first
      do while (associated(property))
         select type (property)
         class is (type_real_property); call field%set_attribute(property%name,property%value)
         end select
         property => property%next
      end do
   end function register

   end subroutine register_fabm_variables
!EOC

!-----------------------------------------------------------------------

   end module register_all_variables

!-----------------------------------------------------------------------
! Copyright by the Bolding & Bruggeman under the GNU Public License    !
!-----------------------------------------------------------------------
#endif

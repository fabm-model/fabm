#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Helper fort the Python interface to FABM. This functionality may move to the FABM core in time.
!
! !INTERFACE:
   module fabm_python_helper
!
! !DESCRIPTION:
! TODO
!
! !USES:
   use fabm, only: type_model, fabm_link_bulk_data, fabm_link_horizontal_data, fabm_link_scalar_data
   use fabm_types

   implicit none
   
   private

   public get_environment_metadata
!EOP
!-----------------------------------------------------------------------

   contains

   subroutine get_environment_metadata(model,environment_names,environment_units)
     type (type_model),                           intent(inout) :: model
     character(len=1024),dimension(:),allocatable,intent(out)   :: environment_names,environment_units

     integer                         :: n
     type (type_link),       pointer :: link
   
      ! Get number of environmental dependencies (light, temperature, etc.)
      n = 0
      link => model%root%first_link
      do while (associated(link))
         if (link%owner) then
            select type (object=>link%target)
               class is (type_bulk_variable)
                  if (allocated(object%alldata).and.object%state_indices%is_empty()) then
                     if (.not.associated(object%alldata(1)%p%p)) n = n+1
                  end if
               class is (type_horizontal_variable)
                  if (allocated(object%alldata).and.object%state_indices%is_empty()) then
                     if (.not.associated(object%alldata(1)%p%p)) n = n+1
                  end if
               class is (type_scalar_variable)
                  if (allocated(object%alldata).and.object%state_indices%is_empty()) then
                     if (.not.associated(object%alldata(1)%p%p)) n = n+1
                  end if
            end select
         end if
         link => link%next
      end do
      
      ! Allocate arrays to hold information on environment
      allocate(environment_names(n))
      allocate(environment_units(n))
      
      ! Get metadata on environmental dependencies (light, temperature, etc.)
      n = 0
      link => model%root%first_link
      do while (associated(link))
         if (link%owner) then
            select type (object=>link%target)
               class is (type_bulk_variable)
                  if (allocated(object%alldata).and.object%state_indices%is_empty()) then
                     if (.not.associated(object%alldata(1)%p%p)) then
                        n = n + 1
                        if (object%standard_variable%is_null()) then
                           environment_names(n) = trim(link%name)
                           environment_units(n) = trim(object%units)
                        else
                           environment_names(n) = trim(object%standard_variable%name)
                           environment_units(n) = trim(object%standard_variable%units)
                        end if
                     end if
                  end if
               class is (type_horizontal_variable)
                  if (allocated(object%alldata).and.object%state_indices%is_empty()) then
                     if (.not.associated(object%alldata(1)%p%p)) then
                        n = n + 1
                        if (object%standard_variable%is_null()) then
                           environment_names(n) = trim(link%name)
                           environment_units(n) = trim(object%units)
                        else
                           environment_names(n) = trim(object%standard_variable%name)
                           environment_units(n) = trim(object%standard_variable%units)
                        end if
                     end if
                  end if
               class is (type_scalar_variable)
                  if (allocated(object%alldata).and.object%state_indices%is_empty()) then
                     if (.not.associated(object%alldata(1)%p%p)) then
                        n = n + 1
                        if (object%standard_variable%is_null()) then
                           environment_names(n) = trim(link%name)
                           environment_units(n) = trim(object%units)
                        else
                           environment_names(n) = trim(object%standard_variable%name)
                           environment_units(n) = trim(object%standard_variable%units)
                        end if
                     end if
                  end if
            end select
         end if
         link => link%next
      end do
   end subroutine

   end module fabm_python_helper

!-----------------------------------------------------------------------
! Copyright by the FABM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

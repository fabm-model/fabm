#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Helper for the Python interface to FABM. This functionality may move to the FABM core in time.
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

   public get_environment_metadata, get_couplings, get_suitable_masters
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
      link => model%links_postcoupling%first
      do while (associated(link))
         if (.not.link%target%read_indices%is_empty().and.link%target%state_indices%is_empty()) then
            select case (link%target%domain)
               case (domain_bulk)
                  if (.not.associated(model%environment%data(link%target%read_indices%pointers(1)%p)%p)) n = n+1
               case (domain_bottom,domain_surface,domain_horizontal)
                  if (.not.associated(model%environment%data_hz(link%target%read_indices%pointers(1)%p)%p)) n = n+1
               case (domain_scalar)
                  if (.not.associated(model%environment%data_scalar(link%target%read_indices%pointers(1)%p)%p)) n = n+1
            end select
         end if
         link => link%next
      end do

      ! Allocate arrays to hold information on environment
      allocate(environment_names(n))
      allocate(environment_units(n))

      ! Get metadata on environmental dependencies (light, temperature, etc.)
      n = 0
      link => model%links_postcoupling%first
      do while (associated(link))
         if (.not.link%target%read_indices%is_empty().and.link%target%state_indices%is_empty()) then
            select case (link%target%domain)
               case (domain_bulk)
                  if (.not.associated(model%environment%data(link%target%read_indices%pointers(1)%p)%p)) then
                     n = n + 1
                     if (.not.associated(link%target%standard_variable)) then
                        environment_names(n) = trim(link%name)
                        environment_units(n) = trim(link%target%units)
                     else
                        environment_names(n) = trim(link%target%standard_variable%name)
                        environment_units(n) = trim(link%target%standard_variable%units)
                     end if
                  end if
               case (domain_bottom,domain_surface,domain_horizontal)
                  if (.not.associated(model%environment%data_hz(link%target%read_indices%pointers(1)%p)%p)) then
                     n = n + 1
                     if (.not.associated(link%target%standard_variable)) then
                        environment_names(n) = trim(link%name)
                        environment_units(n) = trim(link%target%units)
                     else
                        environment_names(n) = trim(link%target%standard_variable%name)
                        environment_units(n) = trim(link%target%standard_variable%units)
                     end if
                  end if
               case (domain_scalar)
                  if (.not.associated(model%environment%data_scalar(link%target%read_indices%pointers(1)%p)%p)) then
                     n = n + 1
                     if (.not.associated(link%target%standard_variable)) then
                        environment_names(n) = trim(link%name)
                        environment_units(n) = trim(link%target%units)
                     else
                        environment_names(n) = trim(link%target%standard_variable%name)
                        environment_units(n) = trim(link%target%standard_variable%units)
                     end if
                  end if
            end select
         end if
         link => link%next
      end do
   end subroutine get_environment_metadata

   subroutine get_couplings(model,link_list)
      type (type_model),    intent(inout) :: model
      type (type_link_list),intent(inout) :: link_list

      type (type_link),pointer :: link,link2

      call link_list%finalize()
      link => model%root%links%first
      do while (associated(link))
         if (link%original%presence/=presence_internal.and..not.link%original%read_indices%is_empty()) then
            link2 => link_list%append(link%target,link%name)
            link2%original => link%original
         end if
         link => link%next
      end do

   end subroutine get_couplings

   function get_suitable_masters(model,slave) result(link_list)
      type (type_model),                     intent(inout) :: model
      type (type_internal_variable), pointer    :: slave
      type (type_link_list),pointer                        :: link_list

      type (type_link),pointer :: link,link2

      allocate(link_list)
      link => model%root%links%first
      do while (associated(link))
         ! Coupled variables cannot serve as master
         if (associated(link%target,link%original) &   ! Uncoupled
             .and..not.associated(link%target,slave) & ! Not self
             .and..not.(link%original%state_indices%is_empty().and..not.slave%state_indices%is_empty()) & ! state variable if slave is state variable
             .and.link%target%domain==slave%domain) then ! And on same domain
            link2 => link_list%append(link%target,link%name)
            link2%original => link%original
         end if
         link => link%next
      end do
   end function get_suitable_masters

   end module fabm_python_helper

!-----------------------------------------------------------------------
! Copyright by the FABM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

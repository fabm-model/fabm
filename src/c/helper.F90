#include "fabm_driver.h"

module fabm_python_helper

   use fabm, only: type_fabm_model
   use fabm_types

   implicit none

   private

   public get_environment_metadata, get_couplings, get_suitable_masters

contains

   subroutine get_environment_metadata(model, list)
     type (type_fabm_model),    intent(inout) :: model
     type (type_variable_list), intent(out)   :: list

     type (type_link), pointer :: link

      ! Get number of environmental dependencies (light, temperature, etc.)
      link => model%links_postcoupling%first
      do while (associated(link))
         if (.not. link%target%read_indices%is_empty() .and. (link%target%source == source_unknown &
            .or. link%target%source == source_external)) then
            select case (link%target%domain)
            case (domain_interior)
               if (.not. associated(model%catalog%interior(link%target%catalog_index)%p)) call list%append(link%target)
            case (domain_bottom, domain_surface, domain_horizontal)
               if (.not. associated(model%catalog%horizontal(link%target%catalog_index)%p)) call list%append(link%target)
            case (domain_scalar)
               if (.not. associated(model%catalog%scalar(link%target%catalog_index)%p)) call list%append(link%target)
            end select
         end if
         link => link%next
      end do
   end subroutine get_environment_metadata

   subroutine get_couplings(model, link_list)
      type (type_fabm_model), intent(inout) :: model
      type (type_link_list),  intent(inout) :: link_list

      type (type_link), pointer :: link, link2

      call link_list%finalize()
      link => model%root%links%first
      do while (associated(link))
         if (link%original%presence /= presence_internal .and. .not. link%original%read_indices%is_empty()) then
            link2 => link_list%append(link%target, link%name)
            link2%original => link%original
         end if
         link => link%next
      end do
   end subroutine get_couplings

   function get_suitable_masters(model, source) result(link_list)
      type (type_fabm_model), intent(inout)  :: model
      type (type_internal_variable), pointer :: source
      type (type_link_list),         pointer :: link_list

      type (type_link), pointer :: link, link2

      allocate(link_list)
      link => model%root%links%first
      do while (associated(link))
         ! Coupled variables cannot serve as master
         if (associated(link%target,link%original) &   ! Uncoupled
             .and. .not. associated(link%target, source) & ! Not self
             .and. .not. (link%original%state_indices%is_empty() .and. .not. source%state_indices%is_empty()) & ! state variable if source is state variable
             .and. link%target%domain == source%domain) then ! And on same domain
            link2 => link_list%append(link%target, link%name)
            link2%original => link%original
         end if
         link => link%next
      end do
   end function get_suitable_masters

end module fabm_python_helper

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

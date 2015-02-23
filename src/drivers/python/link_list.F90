module fabm_c_link_list

   use iso_c_binding, only: c_int, c_ptr, c_f_pointer, c_loc

   use fabm_types

   implicit none

contains

   function link_list_count(plist) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_list_count
      type (c_ptr), intent(in), value :: plist
      integer(c_int)                  :: value

      type (type_link_list),pointer :: list

      call c_f_pointer(plist, list)
      value = list%count()
   end function link_list_count

   function link_list_index(plist,index) bind(c) result(pvariable)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_list_index
      type (c_ptr),  intent(in), value :: plist
      integer(c_int),intent(in), value :: index
      type (c_ptr)                     :: pvariable

      type (type_link_list),pointer :: list
      type (type_link),     pointer :: link
      integer                       :: i

      call c_f_pointer(plist, list)
      link => list%first
      do i=2,index
         link => link%next
      end do
      pvariable = c_loc(link%target)
   end function link_list_index

   subroutine link_list_finalize(plist) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_list_finalize
      type (c_ptr),  intent(in), value :: plist

      type (type_link_list),pointer :: list

      call c_f_pointer(plist, list)
      call list%finalize()
      deallocate(list)
   end subroutine link_list_finalize

end module fabm_c_link_list
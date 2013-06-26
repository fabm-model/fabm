#include "fabm_driver.h"

! This module contains types and procedures for generic properties
! that internally can be of different types (real, integer, logical)
! It also offers a type that represents a property dictionary.

module fabm_properties

   implicit none

   private

#ifdef _FABM_F2003_

   public type_property,type_integer_property,type_real_property,type_logical_property,type_property_dictionary

   integer, parameter :: rk = _FABM_REAL_KIND_
   integer, parameter :: property_string_length = 256

   type type_property
      character(len=property_string_length) :: name       = ''
      character(len=property_string_length) :: long_name  = ''
      character(len=property_string_length) :: units      = ''
      class (type_property), pointer        :: next       => null()
   end type

   type,extends(type_property) :: type_integer_property
      integer :: value
   end type

   type,extends(type_property) :: type_real_property
      real(rk) :: value
   end type

   type,extends(type_property) :: type_logical_property
      logical :: value
   end type

   type type_property_dictionary
      class (type_property),pointer :: first => null()

      contains

      procedure :: set_property
      procedure :: set_real
      procedure :: set_integer
      procedure :: set_logical

      procedure :: get_property
      procedure :: get_real
      procedure :: get_integer
      procedure :: get_logical

      procedure :: update
   end type

contains

   subroutine set_property(dictionary,property,overwrite)
      class (type_property_dictionary),intent(inout) :: dictionary
      class (type_property),           intent(in)    :: property
      logical,optional,                intent(in)    :: overwrite

      class (type_property),pointer                  :: current,next
      logical                                        :: overwrite_eff

      overwrite_eff = .true.
      if (present(overwrite)) overwrite_eff = overwrite

      nullify(next)
      if (.not.associated(dictionary%first)) then
         ! First property in list
         allocate(dictionary%first,source=property)
         current => dictionary%first
      else
         ! List already contains one or more properties.
         if (dictionary%first%name==property%name) then
            ! The provided property replaces the head of the list.
            if (.not.overwrite_eff) return
            next => dictionary%first%next
            deallocate(dictionary%first)
            allocate(dictionary%first,source=property)
            current => dictionary%first
         else
            ! Look for last element in list, or the one prior to the element with the same name.
            current => dictionary%first
            do while (associated(current%next))
               if (current%next%name==property%name) exit
               current => current%next
            end do
            if (associated(current%next)) then
               ! We are replacing current%next
               if (.not.overwrite_eff) return
               next => current%next%next
               deallocate(current%next)
            end if
            allocate(current%next,source=property)
            current => current%next
         end if
      end if
      current%next => next
   end subroutine

   subroutine update(target,source,overwrite)
      class (type_property_dictionary),intent(inout) :: target
      class (type_property_dictionary),intent(in)    :: source
      logical,optional,                intent(in)    :: overwrite
      class (type_property),pointer                  :: current

      current => source%first
      do while (associated(current))
         call target%set_property(current,overwrite)
         current => current%next
      end do
   end subroutine

   subroutine set_real(dictionary,name,value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      real(rk),                        intent(in)    :: value
      call dictionary%set_property(type_real_property(name=name,value=value))
   end subroutine

   subroutine set_integer(dictionary,name,value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      integer,                         intent(in)    :: value
      call dictionary%set_property(type_integer_property(name=name,value=value))
   end subroutine

   subroutine set_logical(dictionary,name,value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      logical,                         intent(in)    :: value
      call dictionary%set_property(type_logical_property(name=name,value=value))
   end subroutine

   function get_property(dictionary,name) result(property)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      class (type_property),pointer                  :: property

      property => dictionary%first
      do while (associated(property))
         if (property%name==name) return
         property => property%next
      end do
   end function

   function get_real(dictionary,name,default) result(value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      real(rk),                        intent(in)    :: default
      real(rk)                                       :: value
      class (type_property),pointer                  :: property
      
      value = default
      property => dictionary%get_property(name)
      if (associated(property)) then
         select type (property)
            class is (type_real_property)
               value = property%value
         end select
      end if
   end function

   function get_integer(dictionary,name,default) result(value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      integer,                         intent(in)    :: default
      integer                                        :: value
      class (type_property),pointer                  :: property
      
      value = default
      property => dictionary%get_property(name)
      if (associated(property)) then
         select type (property)
            class is (type_integer_property)
               value = property%value
         end select
      end if
   end function

   function get_logical(dictionary,name,default) result(value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      logical,                         intent(in)    :: default
      logical                                        :: value
      class (type_property),pointer                  :: property
      
      value = default
      property => dictionary%get_property(name)
      if (associated(property)) then
         select type (property)
            class is (type_logical_property)
               value = property%value
         end select
      end if
   end function
   
#endif

end module fabm_properties
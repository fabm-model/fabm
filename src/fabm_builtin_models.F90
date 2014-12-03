#include "fabm_driver.h"

module fabm_builtin_models
   use fabm_types
   use fabm_standard_variables

   implicit none

   private

   public type_weighted_sum,type_horizontal_weighted_sum,type_simple_depth_integral

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: builtin_factory

   type type_component
      character(len=attribute_length) :: name   = ''
      real(rk)                        :: weight = 1._rk
      logical                         :: include_background = .false.
      type (type_dependency_id)       :: id
      type (type_component),pointer   :: next   => null()
   end type

   type type_horizontal_component
      character(len=attribute_length)          :: name   = ''
      real(rk)                                 :: weight = 1._rk
      logical                                  :: include_background = .false.
      type (type_horizontal_dependency_id)     :: id
      type (type_horizontal_component),pointer :: next   => null()
   end type

   type,extends(type_base_model) :: type_weighted_sum
      character(len=attribute_length) :: output_long_name = ''
      character(len=attribute_length) :: output_units     = ''
      real(rk)                        :: offset           = 0.0_rk
      type (type_diagnostic_variable_id) :: id_output
      type (type_component),pointer   :: first => null()
   contains
      procedure :: initialize     => weighted_sum_initialize
      procedure :: add_component  => weighted_sum_add_component
      procedure :: do             => weighted_sum_do
      procedure :: after_coupling => weighted_sum_after_coupling
      procedure :: add_to_parent  => weighted_sum_add_to_parent
   end type

   type,extends(type_base_model) :: type_horizontal_weighted_sum
      character(len=attribute_length) :: output_long_name = ''
      character(len=attribute_length) :: output_units     = ''
      real(rk)                        :: offset           = 0.0_rk
      type (type_horizontal_diagnostic_variable_id) :: id_output
      type (type_horizontal_component),pointer   :: first => null()
   contains
      procedure :: add_component       => horizontal_weighted_sum_add_component
      procedure :: initialize          => horizontal_weighted_sum_initialize
      procedure :: evaluate_horizontal => horizontal_weighted_sum_evaluate_horizontal
      procedure :: do_bottom           => horizontal_weighted_sum_do_bottom
      procedure :: after_coupling      => horizontal_weighted_sum_after_coupling
      procedure :: add_to_parent       => horizontal_weighted_sum_add_to_parent
   end type

   type,extends(type_base_model) :: type_simple_depth_integral
      type (type_dependency_id)                     :: id_input
      type (type_horizontal_dependency_id)          :: id_depth
      type (type_horizontal_diagnostic_variable_id) :: id_output
      real(rk)                                      :: minimum_depth = 0.0_rk
      real(rk)                                      :: maximum_depth = huge(1.0_rk)
      logical                                       :: average       = .false.
   contains
      procedure :: initialize => simple_depth_integral_initialize
      procedure :: do         => simple_depth_integral_do
   end type

   type,extends(type_base_model) :: type_bulk_constant
      type (type_diagnostic_variable_id) :: id_constant
   contains
      procedure :: initialize => bulk_constant_initialize
   end type

   type,extends(type_base_model) :: type_horizontal_constant
      type (type_horizontal_diagnostic_variable_id) :: id_constant
   contains
      procedure :: initialize => horizontal_constant_initialize
   end type

   contains


   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('bulk_constant');       allocate(type_bulk_constant::model)
         case ('horizontal_constant'); allocate(type_horizontal_constant::model)
         ! Add new examples models here
      end select

   end subroutine

   function weighted_sum_add_to_parent(self,parent,name,create_for_one) result(sum_used)
      class (type_weighted_sum),intent(inout),target :: self
      class (type_base_model),  intent(inout),target :: parent
      character(len=*),         intent(in)           :: name
      logical,optional,         intent(in)           :: create_for_one

      logical :: sum_used,create_for_one_

      create_for_one_ = .false.
      if (present(create_for_one)) create_for_one_ = create_for_one

      call parent%add_bulk_variable(name,self%output_units,name)
      if (.not.associated(self%first)) then
         ! No components - add link to zero field to parent.
         call parent%request_coupling(name,'zero')
         sum_used = .false.
      elseif (.not.associated(self%first%next).and.self%first%weight==1.0_rk.and..not.create_for_one_) then
         ! One component with scale factor 1 - add link to component to parent.
         call parent%request_coupling(name,self%first%name)
         sum_used = .false.
      else
         ! One component with scale factor non equal to 1, or multiple components. Create the sum.
         call parent%add_child(self,trim(name)//'_calculator',configunit=-1)
         call parent%request_coupling(name,trim(name)//'_calculator/result')
         sum_used = .true.
      end if
   end function

   subroutine weighted_sum_initialize(self,configunit)
      class (type_weighted_sum),intent(inout),target :: self
      integer,                  intent(in)           :: configunit

      type (type_component),pointer :: component
      integer           :: i
      character(len=10) :: temp

      i = 0
      component => self%first
      do while (associated(component))
         i = i + 1
         write (temp,'(i0)') i
         call self%register_dependency(component%id,'term'//trim(temp),self%output_units,'term '//trim(temp))
         call self%request_coupling(component%id,trim(component%name))
         component => component%next
      end do
      if (self%output_long_name=='') self%output_long_name = 'result'
      call self%register_diagnostic_variable(self%id_output,'result',self%output_units,'result')
   end subroutine

   subroutine weighted_sum_add_component(self,name,weight,include_background)
      class (type_weighted_sum),intent(inout) :: self
      character(len=*),         intent(in)    :: name
      real(rk),optional,        intent(in)    :: weight
      logical,optional,         intent(in)    :: include_background

      type (type_component),pointer :: component

      if (.not.associated(self%first)) then
         allocate(self%first)
         component => self%first
      else
         component => self%first
         do while (associated(component%next))
            component => component%next
         end do
         allocate(component%next)
         component => component%next
      end if
      component%name = name
      if (present(weight)) component%weight = weight
      if (present(include_background)) component%include_background = include_background
   end subroutine

   subroutine weighted_sum_after_coupling(self)
      class (type_weighted_sum),intent(inout) :: self

      type (type_component),pointer :: component

      ! At this stage, the background values for all variables (if any) are fixed. We can therefore
      ! compute background contributions already, and add those to the space- and time-invariant offset.
      component => self%first
      do while (associated(component))
         if (component%include_background) self%offset = self%offset + component%weight*component%id%background
         component => component%next
      end do
   end subroutine

   subroutine weighted_sum_do(self,_ARGUMENTS_DO_)
      class (type_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      type (type_component),pointer        :: component
      real(rk)                             :: value
      real(rk) _DIMENSION_SLICE_AUTOMATIC_ :: sum

      ! Initialize sum to starting value (typically zero).
      sum = self%offset

      ! Enumerate components included in the sum, and add their contributions.
      component => self%first
      do while (associated(component))
         _LOOP_BEGIN_
            _GET_(component%id,value)
            sum _INDEX_SLICE_ = sum _INDEX_SLICE_ + component%weight*value
         _LOOP_END_
         component => component%next
      end do

      ! Transfer summed values to diagnostic.
      _LOOP_BEGIN_
         _SET_DIAGNOSTIC_(self%id_output,sum _INDEX_SLICE_)
      _LOOP_END_
   end subroutine

   function horizontal_weighted_sum_add_to_parent(self,parent,name,create_for_one) result(sum_used)
      class (type_horizontal_weighted_sum),intent(inout),target :: self
      class (type_base_model),             intent(inout),target :: parent
      character(len=*),                    intent(in)           :: name
      logical,optional,                    intent(in)           :: create_for_one

      logical :: sum_used,create_for_one_

      create_for_one_ = .false.
      if (present(create_for_one)) create_for_one_ = create_for_one

      call parent%add_horizontal_variable(name,self%output_units,name)
      if (.not.associated(self%first)) then
         ! No components - add link to zero field to parent.
         call parent%request_coupling(name,'zero_hz')
         sum_used = .false.
      elseif (.not.associated(self%first%next).and.self%first%weight==1.0_rk.and..not.create_for_one_) then
         ! One component with scale factor 1 - add link to component to parent.
         call parent%request_coupling(name,self%first%name)
         sum_used = .false.
      else
         ! One component with scale factor unequal to 1, or multiple components. Create the sum.
         call parent%add_child(self,trim(name)//'_calculator',configunit=-1)
         call parent%request_coupling(name,trim(name)//'_calculator/result')
         sum_used = .true.
      end if
   end function

   subroutine horizontal_weighted_sum_initialize(self,configunit)
      class (type_horizontal_weighted_sum),intent(inout),target :: self
      integer,                             intent(in)           :: configunit

      type (type_horizontal_component),pointer :: component
      integer           :: i
      character(len=10) :: temp

      i = 0
      component => self%first
      do while (associated(component))
         i = i + 1
         write (temp,'(i0)') i
         call self%register_dependency(component%id,'term'//trim(temp),self%output_units,'term '//trim(temp))
         call self%request_coupling(component%id,trim(component%name))
         component => component%next
      end do
      if (self%output_long_name=='') self%output_long_name = 'result'
      call self%register_diagnostic_variable(self%id_output,'result',self%output_units,'result')
   end subroutine

   subroutine horizontal_weighted_sum_after_coupling(self)
      class (type_horizontal_weighted_sum),intent(inout) :: self

      type (type_horizontal_component),pointer :: component

      ! At this stage, the background values for all variables (if any) are fixed. We can therefore
      ! compute background contributions already, and add those to the space- and time-invariant offset.
      component => self%first
      do while (associated(component))
         if (component%include_background) self%offset = self%offset + component%weight*component%id%background
         component => component%next
      end do
   end subroutine

   subroutine horizontal_weighted_sum_add_component(self,name,weight,include_background)
      class (type_horizontal_weighted_sum),intent(inout) :: self
      character(len=*),                    intent(in)    :: name
      real(rk),optional,                   intent(in)    :: weight
      logical,optional,                    intent(in)    :: include_background

      type (type_horizontal_component),pointer :: component

      if (.not.associated(self%first)) then
         allocate(self%first)
         component => self%first
      else
         component => self%first
         do while (associated(component%next))
            component => component%next
         end do
         allocate(component%next)
         component => component%next
      end if
      component%name = name
      if (present(weight)) component%weight = weight
      if (present(include_background)) component%include_background = include_background
   end subroutine

   subroutine horizontal_weighted_sum_evaluate_horizontal(self,_ARGUMENTS_HZ_)
      class (type_horizontal_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_HZ_

      type (type_horizontal_component),pointer        :: component
      real(rk)                                        :: value
      real(rk) _DIMENSION_HORIZONTAL_SLICE_AUTOMATIC_ :: sum

      sum = self%offset

      component => self%first
      do while (associated(component))
         _HORIZONTAL_LOOP_BEGIN_
            _GET_HORIZONTAL_(component%id,value)
            sum _INDEX_HORIZONTAL_SLICE_ = sum _INDEX_HORIZONTAL_SLICE_ + component%weight*value
         _HORIZONTAL_LOOP_END_
         component => component%next
      end do

      _HORIZONTAL_LOOP_BEGIN_
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output,sum _INDEX_HORIZONTAL_SLICE_)
      _HORIZONTAL_LOOP_END_
   end subroutine

   subroutine horizontal_weighted_sum_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_horizontal_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      call self%evaluate_horizontal(_ARGUMENTS_HZ_)
   end subroutine

   subroutine simple_depth_integral_initialize(self,configunit)
      class (type_simple_depth_integral),intent(inout),target :: self
      integer,                           intent(in)           :: configunit
      call self%register_dependency(self%id_input,'source','','source')
      if (.not.self%average) call self%register_dependency(self%id_depth,standard_variables%bottom_depth)
      call self%register_diagnostic_variable(self%id_output,'result','','result')
   end subroutine

   subroutine simple_depth_integral_do(self,_ARGUMENTS_DO_)
      class (type_simple_depth_integral),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: value,depth

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_input,value)
         if (.not.self%average) then
            _GET_HORIZONTAL_(self%id_depth,depth)
            value = value*(min(self%maximum_depth,depth)-self%minimum_depth)
         end if
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output,value)
      _HORIZONTAL_LOOP_END_
   end subroutine

   subroutine bulk_constant_initialize(self,configunit)
      class (type_bulk_constant),intent(inout),target :: self
      integer,                   intent(in)           :: configunit

      character(len=attribute_length) :: standard_name
      real(rk)                        :: value

      call self%get_parameter(standard_name,'standard_name','','standard name')
      call self%get_parameter(value,'value','','value')
      if (standard_name/='') then
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, standard_variable=type_bulk_standard_variable(name=standard_name))
      else
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none)
      end if
   end subroutine

   subroutine horizontal_constant_initialize(self,configunit)
      class (type_horizontal_constant),intent(inout),target :: self
      integer,                         intent(in)           :: configunit

      character(len=attribute_length) :: standard_name
      real(rk)                        :: value

      call self%get_parameter(standard_name,'standard_name','','standard name')
      call self%get_parameter(value,'value','','value')
      if (standard_name/='') then
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, standard_variable=type_horizontal_standard_variable(name=standard_name))
      else
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none)
      end if
   end subroutine

end module fabm_builtin_models

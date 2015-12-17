#include "fabm_driver.h"

module fabm_builtin_models
   use fabm_types
   use fabm_standard_variables

   implicit none

   private

   public type_weighted_sum,type_horizontal_weighted_sum,type_depth_integral,copy_fluxes

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
      character(len=attribute_length) :: units         = ''
      integer                         :: result_output = output_instantaneous
      real(rk)                        :: offset        = 0.0_rk
      integer                         :: access        = access_read
      type (type_bulk_standard_variable),pointer :: standard_variable => null()
      type (type_diagnostic_variable_id) :: id_output
      type (type_component),pointer   :: first => null()
   contains
      procedure :: initialize     => weighted_sum_initialize
      procedure :: add_component  => weighted_sum_add_component
      procedure :: do             => weighted_sum_do
      procedure :: after_coupling => weighted_sum_after_coupling
      procedure :: add_to_parent  => weighted_sum_add_to_parent
   end type

   type,extends(type_base_model) :: type_weighted_sum_sms_distributor
      type (type_dependency_id)                 :: id_total_sms
      type (type_state_variable_id),allocatable :: id_targets(:)
      real(rk),allocatable                      :: weights(:)
   contains
      !procedure :: do => weighted_sum_sms_distributor_do
   end type

   type,extends(type_base_model) :: type_scaled_interior_variable
      type (type_dependency_id)          :: id_source
      type (type_diagnostic_variable_id) :: id_result
      real(rk)                           :: weight = 1.0_rk
      real(rk)                           :: offset = 0.0_rk
      logical                            :: include_background = .false.
   contains
      procedure :: do             => scaled_interior_variable_do
      procedure :: after_coupling => scaled_interior_variable_after_coupling
   end type

   type,extends(type_base_model) :: type_interior_sms_scaler
      type (type_dependency_id)     :: id_source_sms
      type (type_state_variable_id) :: id_target
      real(rk)                      :: weight = 1.0_rk
   contains
      procedure :: do => interior_sms_scaler_do
   end type

   type,extends(type_base_model) :: type_horizontal_weighted_sum
      character(len=attribute_length) :: units         = ''
      integer                         :: result_output = output_instantaneous
      real(rk)                        :: offset        = 0.0_rk
      integer                         :: access        = access_read
      type (type_horizontal_diagnostic_variable_id) :: id_output
      type (type_horizontal_component),pointer   :: first => null()
   contains
      procedure :: add_component       => horizontal_weighted_sum_add_component
      procedure :: initialize          => horizontal_weighted_sum_initialize
      procedure :: do_horizontal       => horizontal_weighted_sum_do_horizontal
      procedure :: after_coupling      => horizontal_weighted_sum_after_coupling
      procedure :: add_to_parent       => horizontal_weighted_sum_add_to_parent
   end type

   type,extends(type_base_model) :: type_depth_integral
      type (type_dependency_id)                     :: id_input
      type (type_dependency_id)                     :: id_thickness
      type (type_horizontal_diagnostic_variable_id) :: id_output
      real(rk)                                      :: minimum_depth = 0.0_rk
      real(rk)                                      :: maximum_depth = huge(1.0_rk)
      logical                                       :: average       = .false.
   contains
      procedure :: initialize => depth_integral_initialize
      procedure :: get_light  => depth_integral_do_column
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

   type,extends(type_base_model) :: type_constant_surface_flux
      type (type_state_variable_id) :: id_target
      real(rk) :: flux
   contains
      procedure :: initialize => constant_surface_flux_initialize
      procedure :: do_surface => constant_surface_flux_do_surface
   end type

   type,extends(type_base_model) :: type_external_surface_flux
      type (type_state_variable_id)        :: id_target
      type (type_horizontal_dependency_id) :: id_flux
   contains
      procedure :: initialize => external_surface_flux_initialize
      procedure :: do_surface => external_surface_flux_do_surface
   end type

   type,extends(type_base_model) :: type_flux_copier
      type (type_state_variable_id)        :: id_target
      type (type_dependency_id)            :: id_sms
      type (type_horizontal_dependency_id) :: id_bottom_flux, id_surface_flux
   contains
      procedure :: do         => flux_copier_do
      procedure :: do_surface => flux_copier_do_surface
      procedure :: do_bottom  => flux_copier_do_bottom
   end type

   contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('bulk_constant');       allocate(type_bulk_constant::model)
         case ('horizontal_constant'); allocate(type_horizontal_constant::model)
         case ('surface_flux');        allocate(type_constant_surface_flux::model)
         case ('constant_surface_flux');allocate(type_constant_surface_flux::model)
         case ('external_surface_flux');allocate(type_external_surface_flux::model)
         ! Add new examples models here
      end select

   end subroutine

   function weighted_sum_add_to_parent(self,parent,name,create_for_one) result(sum_used)
      class (type_weighted_sum),intent(inout),target :: self
      class (type_base_model),  intent(inout),target :: parent
      character(len=*),         intent(in)           :: name
      logical,optional,         intent(in)           :: create_for_one

      logical                                       :: sum_used,create_for_one_
      type (type_link),                     pointer :: link
      class (type_scaled_interior_variable),pointer :: scaled_variable
      class (type_interior_sms_scaler),     pointer :: sms_scaler

      create_for_one_ = .false.
      if (present(create_for_one)) create_for_one_ = create_for_one

      sum_used = .false.
      if (associated(self%standard_variable)) then
         call parent%add_interior_variable(name,self%units,name,link=link,act_as_state_variable=iand(self%access,access_set_source)/=0,standard_variable=self%standard_variable)
      else
         call parent%add_interior_variable(name,self%units,name,link=link,act_as_state_variable=iand(self%access,access_set_source)/=0)
      end if
      if (.not.associated(self%first)) then
         ! No components - add link to zero field to parent.
         call parent%request_coupling(link,'zero')
      elseif (.not.associated(self%first%next).and.self%first%weight==1.0_rk.and..not.create_for_one_) then
         ! One component with scale factor 1 - add link to component to parent.
         call parent%request_coupling(link,self%first%name)
      elseif (.not.associated(self%first%next)) then
         ! One component with scale factor other than 1 (or a user-specified requirement NOT to make a direct link to the source variable)
         allocate(scaled_variable)
         call parent%add_child(scaled_variable,trim(name)//'_calculator',configunit=-1)
         call scaled_variable%register_dependency(scaled_variable%id_source,'source',self%units,'source variable')
         call scaled_variable%request_coupling(scaled_variable%id_source,self%first%name)
         call scaled_variable%register_diagnostic_variable(scaled_variable%id_result,'result',self%units,'result',output=self%result_output,act_as_state_variable=iand(self%access,access_set_source)/=0)
         scaled_variable%weight = self%first%weight
         scaled_variable%include_background = self%first%include_background
         scaled_variable%offset = self%offset
         call parent%request_coupling(link,trim(name)//'_calculator/result')
         if (iand(self%access,access_set_source)/=0) then
            ! This scaled variable acts as a state variable. Create a child model to distribute source terms to the original source variable.
            allocate(sms_scaler)
            call scaled_variable%add_child(sms_scaler,'sms_distributor',configunit=-1)
            sms_scaler%weight = 1/scaled_variable%weight
            call sms_scaler%register_dependency(sms_scaler%id_source_sms,'source_sms',trim(self%units)//' s-1','sources-sinks of source variable')
            call sms_scaler%request_coupling(sms_scaler%id_source_sms,'result_sms_tot')
            call sms_scaler%register_state_dependency(sms_scaler%id_target,'target',self%units,'target variable')
            call sms_scaler%request_coupling(sms_scaler%id_target,self%first%name)
         end if
      else
         ! Multiple components. Create the sum.
         call parent%add_child(self,trim(name)//'_calculator',configunit=-1)
         call parent%request_coupling(link,trim(name)//'_calculator/result')
         sum_used = .true.
      end if
   end function

   subroutine weighted_sum_initialize(self,configunit)
      class (type_weighted_sum),intent(inout),target :: self
      integer,                  intent(in)           :: configunit

      type (type_component),pointer :: component
      integer           :: ncomponents
      character(len=10) :: temp
      class (type_weighted_sum_sms_distributor), pointer :: sms_distributor

      ncomponents = 0
      component => self%first
      do while (associated(component))
         ncomponents = ncomponents + 1
         write (temp,'(i0)') ncomponents
         call self%register_dependency(component%id,'term'//trim(temp),self%units,'term '//trim(temp))
         call self%request_coupling(component%id,trim(component%name))
         component => component%next
      end do
      call self%register_diagnostic_variable(self%id_output,'result',self%units,'result',output=self%result_output) !,act_as_state_variable=iand(self%access,access_set_source)/=0)

      if (iand(self%access,access_set_source)/=0) then
         ! NB this does not function yet (hence the commented out act_as_state_variable above)
         ! Auto-generation of result_sms_tot fails and the do routine of type_weighted_sum_sms_distributor is not yet implemented.

         ! The sum will act as a state variable. Any source terms will have to be distributed over the individual variables that contribute to the sum.
         allocate(sms_distributor)
         call self%add_child(sms_distributor,'sms_distributor',configunit=-1)
         call sms_distributor%register_dependency(sms_distributor%id_total_sms,'total_sms',trim(self%units)//'/s','sources-sinks of sum')
         call sms_distributor%request_coupling(sms_distributor%id_total_sms,'result_sms_tot')
         allocate(sms_distributor%weights(ncomponents))
         allocate(sms_distributor%id_targets(ncomponents))

         ncomponents = 0
         component => self%first
         do while (associated(component))
            ncomponents = ncomponents + 1
            write (temp,'(i0)') ncomponents
            call sms_distributor%register_state_dependency(sms_distributor%id_targets(ncomponents),'target'//trim(temp),self%units,'target '//trim(temp))
            call sms_distributor%request_coupling(sms_distributor%id_targets(ncomponents),trim(component%name))
            sms_distributor%weights(ncomponents) = component%weight
            component => component%next
         end do
      end if
   end subroutine

   subroutine weighted_sum_add_component(self,name,weight,include_background)
      class (type_weighted_sum),intent(inout) :: self
      character(len=*),         intent(in)    :: name
      real(rk),optional,        intent(in)    :: weight
      logical,optional,         intent(in)    :: include_background

      type (type_component),pointer :: component

      if (_VARIABLE_REGISTERED_(self%id_output)) &
         call self%fatal_error('weighted_sum_add_component','cannot be called after model initialization')

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
      real(rk) :: background

      ! At this stage, the background values for all variables (if any) are fixed. We can therefore
      ! compute background contributions already, and add those to the space- and time-invariant offset.
      background = 0
      component => self%first
      do while (associated(component))
         if (component%include_background) then
            self%offset = self%offset + component%weight*component%id%background
         else
            background = background + component%weight*component%id%background
         end if
         component => component%next
      end do
      call self%id_output%link%target%background_values%set_value(background)
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

   subroutine scaled_interior_variable_do(self,_ARGUMENTS_DO_)
      class (type_scaled_interior_variable),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: value

      _LOOP_BEGIN_
         _GET_(self%id_source,value)
         _SET_DIAGNOSTIC_(self%id_result,self%offset + self%weight*value)
      _LOOP_END_
   end subroutine scaled_interior_variable_do

   subroutine interior_sms_scaler_do(self,_ARGUMENTS_DO_)
      class (type_interior_sms_scaler),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: sms

      _LOOP_BEGIN_
         _GET_(self%id_source_sms,sms)
         _SET_ODE_(self%id_target,sms*self%weight)
      _LOOP_END_
   end subroutine interior_sms_scaler_do

   subroutine scaled_interior_variable_after_coupling(self)
      class (type_scaled_interior_variable),intent(inout) :: self

      if (self%include_background) then
         self%offset = self%offset + self%weight*self%id_source%background
      else
         call self%id_result%link%target%background_values%set_value(self%weight*self%id_source%background)
      end if
   end subroutine scaled_interior_variable_after_coupling

   function horizontal_weighted_sum_add_to_parent(self,parent,name,create_for_one) result(sum_used)
      class (type_horizontal_weighted_sum),intent(inout),target :: self
      class (type_base_model),             intent(inout),target :: parent
      character(len=*),                    intent(in)           :: name
      logical,optional,                    intent(in)           :: create_for_one

      logical :: sum_used,create_for_one_
      type (type_link),pointer :: link

      create_for_one_ = .false.
      if (present(create_for_one)) create_for_one_ = create_for_one

      call parent%add_horizontal_variable(name,self%units,name,link=link,act_as_state_variable=iand(self%access,access_set_source)/=0)
      if (.not.associated(self%first)) then
         ! No components - add link to zero field to parent.
         call parent%request_coupling(link,'zero_hz')
         sum_used = .false.
      elseif (.not.associated(self%first%next).and.self%first%weight==1.0_rk.and..not.create_for_one_) then
         ! One component with scale factor 1 - add link to component to parent.
         call parent%request_coupling(link,self%first%name)
         sum_used = .false.
      else
         ! One component with scale factor unequal to 1, or multiple components. Create the sum.
         call parent%add_child(self,trim(name)//'_calculator',configunit=-1)
         call parent%request_coupling(link,trim(name)//'_calculator/result')
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
         call self%register_dependency(component%id,'term'//trim(temp),self%units,'term '//trim(temp))
         call self%request_coupling(component%id,trim(component%name))
         component => component%next
      end do
      call self%register_diagnostic_variable(self%id_output,'result',self%units,'result',output=self%result_output,source=source_do_horizontal)
   end subroutine

   subroutine horizontal_weighted_sum_after_coupling(self)
      class (type_horizontal_weighted_sum),intent(inout) :: self

      type (type_horizontal_component),pointer :: component
      real(rk) :: background

      ! At this stage, the background values for all variables (if any) are fixed. We can therefore
      ! compute background contributions already, and add those to the space- and time-invariant offset.
      background = 0
      component => self%first
      do while (associated(component))
         if (component%include_background) then
            self%offset = self%offset + component%weight*component%id%background
         else
            background = background + component%weight*component%id%background
         end if
         component => component%next
      end do
      call self%id_output%link%target%background_values%set_value(background)
   end subroutine

   subroutine horizontal_weighted_sum_add_component(self,name,weight,include_background)
      class (type_horizontal_weighted_sum),intent(inout) :: self
      character(len=*),                    intent(in)    :: name
      real(rk),optional,                   intent(in)    :: weight
      logical,optional,                    intent(in)    :: include_background

      type (type_horizontal_component),pointer :: component

      if (_VARIABLE_REGISTERED_(self%id_output)) &
         call self%fatal_error('weighted_sum_add_component','cannot be called after model initialization')

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

   subroutine horizontal_weighted_sum_do_horizontal(self,_ARGUMENTS_HORIZONTAL_)
      class (type_horizontal_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_

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
   end subroutine horizontal_weighted_sum_do_horizontal

   subroutine depth_integral_initialize(self,configunit)
      class (type_depth_integral),intent(inout),target :: self
      integer,                           intent(in)           :: configunit
      call self%register_dependency(self%id_input,'source','','source')
      call self%register_dependency(self%id_thickness,standard_variables%cell_thickness)
      call self%register_diagnostic_variable(self%id_output,'result','','result',source=source_do_column)
   end subroutine depth_integral_initialize

   subroutine depth_integral_do_column(self,_ARGUMENTS_VERTICAL_)
      class (type_depth_integral),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: h,value,cum,depth
      logical :: started

      cum = 0
      depth = 0
      started = self%minimum_depth<=0
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_thickness,h)
         _GET_(self%id_input,value)

         depth = depth + h
         if (.not.started) then
            ! Not yet at minimum depth before
            if (depth>=self%minimum_depth) then
               ! Now crossing minimum depth interface
               started = .true.
               h = depth-self%minimum_depth
            end if
         elseif (depth>self%maximum_depth) then
            ! Now crossing maximum depth interface; subtract part of layer height that is not included
            h = h - (depth-self%maximum_depth)
            _VERTICAL_LOOP_EXIT_
         end if
         cum = cum + h*value
      _VERTICAL_LOOP_END_

      if (.not.self%average) then
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output,cum)
      elseif (depth>self%minimum_depth) then
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output,cum/(min(self%maximum_depth,depth)-self%minimum_depth))
      endif
   end subroutine depth_integral_do_column

   subroutine bulk_constant_initialize(self,configunit)
      class (type_bulk_constant),intent(inout),target :: self
      integer,                   intent(in)           :: configunit

      character(len=attribute_length) :: standard_name
      real(rk)                        :: value

      call self%get_parameter(standard_name,'standard_name','','standard name',default='')
      call self%get_parameter(value,'value','','value')
      if (standard_name/='') then
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, standard_variable=type_bulk_standard_variable(name=standard_name), source=source_none)
      else
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, source=source_none)
      end if
   end subroutine bulk_constant_initialize

   subroutine horizontal_constant_initialize(self,configunit)
      class (type_horizontal_constant),intent(inout),target :: self
      integer,                         intent(in)           :: configunit

      character(len=attribute_length) :: standard_name
      real(rk)                        :: value

      call self%get_parameter(standard_name,'standard_name','','standard name',default='')
      call self%get_parameter(value,'value','','value')
      if (standard_name/='') then
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, standard_variable=type_horizontal_standard_variable(name=standard_name), source=source_none)
      else
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, source=source_none)
      end if
   end subroutine horizontal_constant_initialize

   subroutine constant_surface_flux_initialize(self,configunit)
      class (type_constant_surface_flux),intent(inout),target :: self
      integer,                  intent(in)           :: configunit

      call self%register_state_dependency(self%id_target,'target','UNITS m-3','target variable')
      call self%get_parameter(self%flux,'flux','UNITS m-2 s-1','flux (positive for into water)')
   end subroutine constant_surface_flux_initialize

   subroutine constant_surface_flux_do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_constant_surface_flux), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      _HORIZONTAL_LOOP_BEGIN_
         _SET_SURFACE_EXCHANGE_(self%id_target,self%flux)
      _HORIZONTAL_LOOP_END_
   end subroutine constant_surface_flux_do_surface

   subroutine external_surface_flux_initialize(self,configunit)
      class (type_external_surface_flux),intent(inout),target :: self
      integer,                  intent(in)           :: configunit

      call self%register_state_dependency(self%id_target,'target','UNITS m-3','target variable')
      call self%register_dependency(self%id_flux,'flux','UNITS m-2 s-1','surface flux')
   end subroutine external_surface_flux_initialize

   subroutine external_surface_flux_do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_external_surface_flux), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: flux

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_flux,flux)
         _SET_SURFACE_EXCHANGE_(self%id_target,flux)
      _HORIZONTAL_LOOP_END_
   end subroutine external_surface_flux_do_surface

   subroutine copy_fluxes(source_model,source_variable,target_variable)
      class (type_base_model),           intent(inout), target :: source_model
      type (type_diagnostic_variable_id),intent(in)            :: source_variable
      type (type_state_variable_id),     intent(in)            :: target_variable
      class (type_flux_copier),pointer :: copier

      allocate(copier)
      call source_model%add_child(copier,'redirect_'//trim(target_variable%link%name)//'_fluxes',configunit=-1)
      call copier%register_state_dependency(copier%id_target,'target','','target variable')
      call copier%register_dependency(copier%id_sms,'sms','','sources minus sinks')
      call copier%register_dependency(copier%id_bottom_flux,'bottom_flux','','bottom flux')
      call copier%register_dependency(copier%id_surface_flux,'surface_flux','','surface flux')
      call copier%request_coupling(copier%id_target,target_variable%link%target%name)
      call copier%request_coupling(copier%id_sms,trim(source_variable%link%target%name)//'_sms_tot')
      call copier%request_coupling(copier%id_bottom_flux,trim(source_variable%link%target%name)//'_bfl_tot')
      call copier%request_coupling(copier%id_surface_flux,trim(source_variable%link%target%name)//'_sfl_tot')
   end subroutine

   subroutine flux_copier_do(self,_ARGUMENTS_DO_)
      class (type_flux_copier), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: sms

      _LOOP_BEGIN_
         _GET_(self%id_sms,sms)
         _SET_ODE_(self%id_target,sms)
      _LOOP_END_
   end subroutine flux_copier_do

   subroutine flux_copier_do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_flux_copier), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: flux

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_surface_flux,flux)
         _SET_SURFACE_EXCHANGE_(self%id_target,flux)
      _HORIZONTAL_LOOP_END_
   end subroutine flux_copier_do_surface

   subroutine flux_copier_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_flux_copier), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: flux

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_bottom_flux,flux)
         _SET_BOTTOM_EXCHANGE_(self%id_target,flux)
      _HORIZONTAL_LOOP_END_
   end subroutine flux_copier_do_bottom

end module fabm_builtin_models

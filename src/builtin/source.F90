#include "fabm_driver.h"

module fabm_builtin_source
   use fabm_types

   implicit none

   private

   public type_interior_source, type_bottom_source, type_surface_source
   public type_constant_surface_flux, type_external_surface_flux, type_external_bottom_flux
   public copy_fluxes, copy_horizontal_fluxes

   type, extends(type_base_model) :: type_constant_surface_flux
      type (type_state_variable_id) :: id_target
      real(rk) :: flux
   contains
      procedure :: initialize => constant_surface_flux_initialize
      procedure :: do_surface => constant_surface_flux_do_surface
   end type

   type, extends(type_base_model) :: type_horizontal_flux
      type (type_link), pointer            :: target => null()
      type (type_horizontal_add_id)        :: id_target_flux
      type (type_horizontal_dependency_id) :: id_flux
      real(rk)                             :: scale_factor = 1.0_rk
   contains
      procedure :: do_horizontal => horizontal_flux_do_horizontal
   end type

   type, extends(type_horizontal_flux) :: type_external_surface_flux
   contains
      procedure :: initialize => external_surface_flux_initialize
   end type

   type, extends(type_horizontal_flux) :: type_external_bottom_flux
   contains
      procedure :: initialize => external_bottom_flux_initialize
   end type

   type, extends(type_base_model) :: type_interior_source
      type (type_link), pointer :: target => null()
      type (type_add_id)        :: id_target_sms
      type (type_dependency_id) :: id_source
      real(rk)                  :: scale_factor = 1.0_rk
   contains
      procedure :: initialize => interior_source_initialize
      procedure :: do         => interior_source_do
   end type

   type, extends(type_horizontal_flux) :: type_bottom_source
   contains
      procedure :: initialize => bottom_source_initialize
   end type

   type, extends(type_horizontal_flux) :: type_surface_source
   contains
      procedure :: initialize => surface_source_initialize
   end type

   interface copy_fluxes
      module procedure copy_fluxes_to_id
      module procedure copy_fluxes_to_named_variable
   end interface

contains

   subroutine constant_surface_flux_initialize(self, configunit)
      class (type_constant_surface_flux), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      call self%register_state_dependency(self%id_target, 'target', 'UNITS m-3', 'target variable')
      call self%get_parameter(self%flux, 'flux', 'UNITS m-2 s-1', 'flux (positive for into water)')
   end subroutine constant_surface_flux_initialize

   subroutine constant_surface_flux_do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_constant_surface_flux), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      _HORIZONTAL_LOOP_BEGIN_
         _ADD_SURFACE_FLUX_(self%id_target, self%flux)
      _HORIZONTAL_LOOP_END_
   end subroutine constant_surface_flux_do_surface

   subroutine horizontal_flux_do_horizontal(self, _ARGUMENTS_HORIZONTAL_)
      class (type_horizontal_flux), intent(in) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_

      real(rk) :: flux

      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_flux, flux)
         _ADD_HORIZONTAL_(self%id_target_flux, self%scale_factor * flux)
      _HORIZONTAL_LOOP_END_
   end subroutine horizontal_flux_do_horizontal

   subroutine external_surface_flux_initialize(self, configunit)
      class (type_external_surface_flux), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      call self%add_interior_variable('target', 'UNITS m-3', 'target variable', link=self%target)
      call self%register_surface_flux(self%target, self%id_target_flux, source=source_do_horizontal)
      call self%register_dependency(self%id_flux, 'flux', 'UNITS m-2 s-1', 'surface flux')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine external_surface_flux_initialize

   subroutine external_bottom_flux_initialize(self, configunit)
      class (type_external_bottom_flux), intent(inout), target :: self
      integer,                           intent(in)            :: configunit

      call self%add_interior_variable('target', 'UNITS m-3', 'target variable', link=self%target)
      call self%register_bottom_flux(self%target, self%id_target_flux, source=source_do_horizontal)
      call self%register_dependency(self%id_flux, 'flux', 'UNITS m-2 s-1', 'bottom flux')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine external_bottom_flux_initialize

   subroutine interior_source_initialize(self, configunit)
      class (type_interior_source), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      call self%add_interior_variable('target', 'UNITS m-3', 'target variable', link=self%target)
      call self%register_source(self%target, self%id_target_sms)
      call self%register_dependency(self%id_source, 'source', 'UNITS m-3 s-1', 'source')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine interior_source_initialize

   subroutine interior_source_do(self, _ARGUMENTS_DO_)
      class (type_interior_source), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: source

      _LOOP_BEGIN_
         _GET_(self%id_source,source)
         _ADD_(self%id_target_sms, self%scale_factor * source)
      _LOOP_END_
   end subroutine interior_source_do

   subroutine bottom_source_initialize(self, configunit)
      class (type_bottom_source), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      call self%add_horizontal_variable('target', 'UNITS m-2', 'target variable', domain=domain_bottom, link=self%target)
      call self%register_bottom_source(self%target, self%id_target_flux, source=source_do_horizontal)
      call self%register_dependency(self%id_flux, 'source', 'UNITS m-2 s-1', 'source')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine bottom_source_initialize

   subroutine surface_source_initialize(self, configunit)
      class (type_surface_source), intent(inout), target :: self
      integer,                     intent(in)            :: configunit

      call self%add_horizontal_variable('target', 'UNITS m-2', 'target variable', domain=domain_bottom, link=self%target)
      call self%register_surface_source(self%target, self%id_target_flux, source=source_do_horizontal)
      call self%register_dependency(self%id_flux, 'source', 'UNITS m-2 s-1', 'source')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine surface_source_initialize

   subroutine copy_fluxes_to_id(source_model, source_variable, target_variable, scale_factor)
      class (type_base_model),            intent(inout), target :: source_model
      type (type_diagnostic_variable_id), intent(in)            :: source_variable
      type (type_state_variable_id),      intent(in)            :: target_variable
      real(rk), optional,                 intent(in)            :: scale_factor
      call copy_fluxes_to_named_variable(source_model, source_variable, target_variable%link%target%name, scale_factor)
   end subroutine

   subroutine copy_fluxes_to_named_variable(source_model, source_variable, target_variable, scale_factor)
      class (type_base_model),            intent(inout), target :: source_model
      type (type_diagnostic_variable_id), intent(in)            :: source_variable
      character(len=*),                   intent(in)            :: target_variable
      real(rk), optional,                 intent(in)            :: scale_factor

      class (type_interior_source),       pointer :: copier
      class (type_external_bottom_flux),  pointer :: bfl_copier
      class (type_external_surface_flux), pointer :: sfl_copier

      allocate(copier)
      call source_model%add_child(copier, '*')
      if (present(scale_factor)) copier%scale_factor = scale_factor
      call copier%request_coupling(copier%target, target_variable)
      call copier%request_coupling(copier%id_source, trim(source_variable%link%target%name) // '_sms_tot')

      allocate(bfl_copier)
      call source_model%add_child(bfl_copier, '*')
      if (present(scale_factor)) bfl_copier%scale_factor = scale_factor
      call bfl_copier%request_coupling(bfl_copier%target, target_variable)
      call bfl_copier%request_coupling(bfl_copier%id_flux, trim(source_variable%link%target%name) // '_bfl_tot')

      allocate(sfl_copier)
      call source_model%add_child(sfl_copier, '*')
      if (present(scale_factor)) sfl_copier%scale_factor = scale_factor
      call sfl_copier%request_coupling(sfl_copier%target, target_variable)
      call sfl_copier%request_coupling(sfl_copier%id_flux, trim(source_variable%link%target%name) // '_sfl_tot')
   end subroutine

   subroutine copy_horizontal_fluxes(source_model, source_variable, target_variable, scale_factor)
      class (type_base_model),                       intent(inout), target :: source_model
      type (type_horizontal_diagnostic_variable_id), intent(in)            :: source_variable
      character(len=*),                              intent(in)            :: target_variable
      real(rk),optional,                             intent(in)            :: scale_factor

      class (type_bottom_source),  pointer :: bottom_copier
      class (type_surface_source), pointer :: surface_copier

      select case (source_variable%link%target%domain)
      case (domain_bottom)
         allocate(bottom_copier)
         call source_model%add_child(bottom_copier, '*')
         if (present(scale_factor)) bottom_copier%scale_factor = scale_factor
         call bottom_copier%request_coupling(bottom_copier%target, target_variable)
         call bottom_copier%request_coupling(bottom_copier%id_flux, trim(source_variable%link%target%name) // '_sms_tot')
      case (domain_surface)
         allocate(surface_copier)
         call source_model%add_child(surface_copier, '*')
         if (present(scale_factor)) surface_copier%scale_factor = scale_factor
         call surface_copier%request_coupling(surface_copier%target, target_variable)
         call surface_copier%request_coupling(surface_copier%id_flux, trim(source_variable%link%target%name) // '_sms_tot')
      case default
         call source_model%fatal_error('copy_horizontal_fluxes','source variable has unknown domain (should be either surface or bottom)')
      end select
   end subroutine copy_horizontal_fluxes

end module

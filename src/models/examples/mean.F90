#include "fabm_driver.h"

module examples_mean

   use fabm_types
   use fabm_expressions

   implicit none

   private

   type,extends(type_base_model),public :: type_examples_mean
      type (type_dependency_id)                  :: id_temp
      type (type_dependency_id)                  :: id_temp_tempmean
      type (type_horizontal_dependency_id)       :: id_temp_vertmean,id_temp_vertmean_20m,id_temp_vertmean_tempmean
      type (type_diagnostic_variable_id)         :: id_temp_tempmean_diag
      type (type_surface_diagnostic_variable_id) :: id_temp_vertmean_diag,id_temp_vertmean_20m_diag,id_temp_vertmean_tempmean_diag
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_examples_mean), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_temp_vertmean, vertical_mean(self%id_temp))
      call self%register_dependency(self%id_temp_vertmean_20m, vertical_mean(self%id_temp, maximum_depth=20._rk))
      call self%register_dependency(self%id_temp_tempmean, temporal_mean(self%id_temp, period=14._rk*86400._rk, resolution=3600._rk))
      call self%register_dependency(self%id_temp_vertmean_tempmean, temporal_mean(self%id_temp_vertmean, period=30._rk*86400._rk, resolution=86400._rk))

      call self%register_diagnostic_variable(self%id_temp_tempmean_diag, 'temp_14dmean',  'degree_C', '14-day running mean of temperature')
      call self%register_diagnostic_variable(self%id_temp_vertmean_diag, 'temp_vertmean', 'degree_C', 'vertical mean temperature')
      call self%register_diagnostic_variable(self%id_temp_vertmean_20m_diag, 'temp_vertmean_20m', 'degree_C', 'vertical mean temperature above 20 m')
      call self%register_diagnostic_variable(self%id_temp_vertmean_tempmean_diag, 'temp_vertmean_tempmean', 'degree_C','30-day running mean of vertical mean temperature')
   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_examples_mean), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: temp

      _LOOP_BEGIN_
         _GET_(self%id_temp_tempmean, temp)
         _SET_DIAGNOSTIC_(self%id_temp_tempmean_diag ,temp)
      _LOOP_END_
   end subroutine do

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_examples_mean), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: temp, temp_20m, temp_tempmean

      _SURFACE_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_temp_vertmean, temp)
         _GET_HORIZONTAL_(self%id_temp_vertmean_20m, temp_20m)
         _GET_HORIZONTAL_(self%id_temp_vertmean_tempmean, temp_tempmean)
         _SET_SURFACE_DIAGNOSTIC_(self%id_temp_vertmean_diag, temp)
         _SET_SURFACE_DIAGNOSTIC_(self%id_temp_vertmean_20m_diag, temp_20m)
         _SET_SURFACE_DIAGNOSTIC_(self%id_temp_vertmean_tempmean_diag, temp_tempmean)
      _SURFACE_LOOP_END_
   end subroutine do_surface

end module

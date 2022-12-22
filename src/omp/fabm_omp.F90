#include "fabm_driver.h"
#include "fabm_private.h"

module fabm_omp

   use omp_lib, only: omp_get_max_threads, omp_get_thread_num

   use fabm
   use fabm_types
   use fabm_work
   use fabm_job, only: type_job, type_task
   use fabm_properties, only: type_property_dictionary

   implicit none

   public

   type type_generic_cache
      type (type_interior_cache)   :: interior
      type (type_horizontal_cache) :: horizontal
      type (type_vertical_cache)   :: vertical
   end type

   type, extends(type_fabm_model) :: type_fabm_omp_model
      type (type_generic_cache), allocatable :: caches(:)
   contains
      procedure :: start
      procedure :: create_cache
      procedure :: process_job
      procedure :: get_interior_cache
      procedure :: get_horizontal_cache
      procedure :: get_vertical_cache
   end type

contains

   function fabm_create_omp_model(path, initialize, parameters, unit) result(model)
      use fabm_config, only: fabm_configure_model
      character(len=*),                optional, intent(in) :: path
      logical,                         optional, intent(in) :: initialize
      type (type_property_dictionary), optional, intent(in) :: parameters
      integer,                         optional, intent(in) :: unit
      class (type_fabm_omp_model), pointer                  :: model

      logical :: initialize_

      ! Make sure the library is initialized.
      call fabm_initialize_library()

      allocate(model)
      call fabm_configure_model(model%root, model%schedules, model%log, path, parameters=parameters, unit=unit)

      ! Initialize model tree
      initialize_ = .true.
      if (present(initialize)) initialize_ = initialize
      if (initialize_) call model%initialize()
   end function fabm_create_omp_model

   subroutine create_cache(self, cache)
      class (type_fabm_omp_model), intent(in)    :: self
      type (type_generic_cache),   intent(inout) :: cache

      call cache_create(self%domain, self%cache_fill_values, cache%interior)
      call cache_create(self%domain, self%cache_fill_values, cache%horizontal)
      call cache_create(self%domain, self%cache_fill_values, cache%vertical)
   end subroutine

   subroutine start(self)
      class (type_fabm_omp_model), intent(inout), target :: self

      integer :: nthread, ithread

      call self%type_fabm_model%start()

      nthread = omp_get_max_threads()
      allocate(self%caches(nthread))
      do ithread = 1, nthread
         call self%create_cache(self%caches(ithread))
      end do
   end subroutine

   function get_interior_cache(self) result(cache)
      class (type_fabm_omp_model), intent(in), target :: self
      type (type_interior_cache), pointer :: cache
      cache => self%caches(omp_get_thread_num() + 1)%interior
   end function

   function get_horizontal_cache(self) result(cache)
      class (type_fabm_omp_model), intent(in), target :: self
      type (type_horizontal_cache), pointer :: cache
      cache => self%caches(omp_get_thread_num() + 1)%horizontal
   end function

   function get_vertical_cache(self) result(cache)
      class (type_fabm_omp_model), intent(in), target :: self
      type (type_vertical_cache), pointer :: cache
      cache => self%caches(omp_get_thread_num() + 1)%vertical
   end function

   subroutine process_job(self, job _POSTARG_HORIZONTAL_LOCATION_RANGE_)
      class (type_fabm_omp_model), intent(inout), target :: self
      type (type_job),             intent(in)            :: job
      _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_RANGE_

      integer                   :: i
      type (type_task), pointer :: task
      _DECLARE_LOCATION_
      integer                   :: ithread

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      ! Jobs must be applied across the entire depth range (if any),
      ! so we set vertical start and stop indices here.
      integer :: _VERTICAL_START_, _VERTICAL_STOP_
      _VERTICAL_START_ = self%domain%start(_FABM_DEPTH_DIMENSION_INDEX_)
      _VERTICAL_STOP_ = self%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_)
#endif

      do i = 1, size(job%interior_store_prefill)
         if (job%interior_store_prefill(i)) then
            _BEGIN_OUTER_INTERIOR_LOOP_
               self%store%interior _INDEX_GLOBAL_INTERIOR_PLUS_1_(_START_:_STOP_, i) = self%store%interior_fill_value(i)
            _END_OUTER_INTERIOR_LOOP_
         end if
      end do
      do i = 1, size(job%horizontal_store_prefill)
         if (job%horizontal_store_prefill(i)) then
            _BEGIN_OUTER_HORIZONTAL_LOOP_
               self%store%horizontal _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_START_:_STOP_, i) = self%store%horizontal_fill_value(i)
            _END_OUTER_HORIZONTAL_LOOP_
         end if
      end do

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ithread,task _POSTARG_LOCATION_)
      ithread = omp_get_thread_num() + 1

      task => job%first_task
      do while (associated(task))
         select case (task%operation)
         case (source_do)
#if (_FABM_DIMENSION_COUNT_==0||(_FABM_DIMENSION_COUNT_==1&&_FABM_VECTORIZED_DIMENSION_INDEX_==1))
            !$OMP SINGLE
#else
            !$OMP DO
#endif
            _BEGIN_OUTER_INTERIOR_LOOP_
               call process_interior_slice(task,self%domain, self%catalog, self%cache_fill_values, self%store, self%caches(ithread)%interior _POSTARG_INTERIOR_IN_)
            _END_OUTER_INTERIOR_LOOP_
#if (_FABM_DIMENSION_COUNT_==0||(_FABM_DIMENSION_COUNT_==1&&_FABM_VECTORIZED_DIMENSION_INDEX_==1))
            !$OMP END SINGLE
#else
            !$OMP END DO
#endif


         case (source_do_surface, source_do_bottom, source_do_horizontal)
#if (_HORIZONTAL_DIMENSION_COUNT_==0||(_HORIZONTAL_DIMENSION_COUNT_==1&&_FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_))
            !$OMP SINGLE
#else
            !$OMP DO
#endif
            _BEGIN_OUTER_HORIZONTAL_LOOP_
               call process_horizontal_slice(task, self%domain, self%catalog, self%cache_fill_values, self%store, self%caches(ithread)%horizontal _POSTARG_HORIZONTAL_IN_)
            _END_OUTER_HORIZONTAL_LOOP_
#if (_HORIZONTAL_DIMENSION_COUNT_==0||(_HORIZONTAL_DIMENSION_COUNT_==1&&_FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_))
            !$OMP END SINGLE
#else
            !$OMP END DO
#endif

         case (source_do_column)
#if (_HORIZONTAL_DIMENSION_COUNT_==0)
            !$OMP SINGLE
#else
#  if _FABM_BOTTOM_INDEX_==-1
#    ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
            !$OMP DO PRIVATE(_VERTICAL_START_)
#    else
            !$OMP DO PRIVATE(_VERTICAL_STOP_)
#    endif
#  else
            !$OMP DO
#  endif
#endif
            _BEGIN_OUTER_VERTICAL_LOOP_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==-1
#    ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
               _VERTICAL_START_ = self%domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#    else
               _VERTICAL_STOP_ = self%domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#    endif
#  endif
#endif
               if (_IS_UNMASKED_(self%domain%mask_hz _INDEX_HORIZONTAL_LOCATION_)) call process_vertical_slice(task, self%domain, &
                  self%catalog, self%cache_fill_values, self%store, self%caches(ithread)%vertical _POSTARG_VERTICAL_IN_)
            _END_OUTER_VERTICAL_LOOP_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
            _VERTICAL_START_ = self%domain%start(_FABM_DEPTH_DIMENSION_INDEX_)
            _VERTICAL_STOP_ = self%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_)
#endif
#if (_HORIZONTAL_DIMENSION_COUNT_==0)
            !$OMP END SINGLE
#endif
         end select
         task => task%next
      end do

      !$OMP END PARALLEL

   end subroutine process_job

end module fabm_omp

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - Public License - www.gnu.org
!-----------------------------------------------------------------------

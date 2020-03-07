#include "fabm_driver.h"
#include "fabm_private.h"

module fabm_work

   use fabm_types, only: rki => rk, rke, type_interior_cache, type_horizontal_cache, type_vertical_cache, &
      prefill_none, prefill_constant, source_do, source_do_horizontal, source_do_surface, source_do_bottom, &
      source_get_light_extinction, source_get_vertical_movement, source2string, &
      domain_interior, domain_surface, domain_bottom, domain_horizontal
   use fabm_job, only: type_task, type_call
   use fabm_driver, only: driver

   implicit none

   private

   public type_domain, type_catalog, type_store, type_cache_fill_values
   public cache_create, cache_pack, cache_unpack
   public process_interior_slice, process_horizontal_slice, process_vertical_slice

   type type_domain
      ! Information about the model domain
      integer :: size(_FABM_DIMENSION_COUNT_)
      integer :: horizontal_size(_HORIZONTAL_DIMENSION_COUNT_)

#ifdef _HAS_MASK_
#  ifndef _FABM_HORIZONTAL_MASK_
      _FABM_MASK_TYPE_, pointer _DIMENSION_GLOBAL_ :: mask => null()
#  endif
      _FABM_MASK_TYPE_, pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: mask_hz => null()
#endif

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==0
      integer :: bottom_index = -1
#  elif _FABM_BOTTOM_INDEX_==-1
      integer, pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: bottom_indices => null()
#  endif
      integer :: surface_index = -1
#endif
   end type

   ! --------------------------------------------------------------------------
   ! Derived types for catalog with pointers to all available fields
   ! --------------------------------------------------------------------------

   type type_interior_data_pointer
      real(rke), pointer _DIMENSION_GLOBAL_ :: p => null()
   end type

   type type_horizontal_data_pointer
      real(rke), pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: p => null()
   end type

   type type_scalar_data_pointer
      real(rke), pointer :: p => null()
   end type

   type type_catalog
      type (type_interior_data_pointer),   allocatable :: interior(:)
      type (type_horizontal_data_pointer), allocatable :: horizontal(:)
      type (type_scalar_data_pointer),     allocatable :: scalar(:)
      integer, allocatable :: interior_sources(:)
      integer, allocatable :: horizontal_sources(:)
      integer, allocatable :: scalar_sources(:)
   end type

   ! --------------------------------------------------------------------------
   ! Derived types for variable store
   ! (spatially explicit model outputs needed by other BGC modules or host)
   ! --------------------------------------------------------------------------

   type type_store
      real(rke), allocatable _DIMENSION_GLOBAL_PLUS_1_            :: interior
      real(rke), allocatable _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: horizontal
      real(rke), allocatable                                      :: interior_fill_value(:)
      real(rke), allocatable                                      :: horizontal_fill_value(:)
      real(rke), allocatable                                      :: interior_missing_value(:)
      real(rke), allocatable                                      :: horizontal_missing_value(:)
   end type

   type type_cache_fill_values
      real(rki), allocatable :: read(:)
      real(rki), allocatable :: read_hz(:)
      real(rki), allocatable :: read_scalar(:)
      real(rki), allocatable :: write(:)
      real(rki), allocatable :: write_hz(:)
   end type

   interface allocate_and_fill
      module procedure allocate_and_fill_0d
      module procedure allocate_and_fill_1d
   end interface

   interface cache_create
      module procedure create_interior_cache
      module procedure create_horizontal_cache
      module procedure create_vertical_cache
   end interface

   interface cache_pack
      module procedure begin_interior_task
      module procedure begin_horizontal_task
      module procedure begin_vertical_task
   end interface

   interface cache_unpack
      module procedure end_interior_task
      module procedure end_horizontal_task
      module procedure end_vertical_task
   end interface

   integer, parameter :: array_block_size = 8

   real(rki), parameter :: not_written = huge(1.0_rki)

contains

   subroutine allocate_and_fill_0d(target, fill, lb, n)
      real(rki), allocatable, intent(out) :: target(:)
      real(rki),              intent(in)  :: fill(:)
      integer,                intent(in)  :: lb, n

      allocate(target(lb:size(fill)))
      target(1:) = fill
   end subroutine

   subroutine allocate_and_fill_1d(target, fill, lb, n)
      real(rki), allocatable, intent(out) :: target(:, :)
      real(rki),              intent(in)  :: fill(:)
      integer,                intent(in)  :: lb, n

      integer :: i

      allocate(target(1:n, lb:size(fill)))
      do i = 1, size(fill)
         target(:, i) = fill(i)
      end do
   end subroutine

   subroutine create_interior_cache(domain, fill_values, cache)
      type (type_domain),            intent(in)  :: domain
      type (type_cache_fill_values), intent(in)  :: fill_values
      type (type_interior_cache),    intent(out) :: cache

      integer :: n, n_mod, i

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      n = domain%size(_FABM_VECTORIZED_DIMENSION_INDEX_)
      n_mod = mod(n, array_block_size)
      if (n_mod /= 0) n = n - n_mod + array_block_size
#  ifdef _HAS_MASK_
      allocate(cache%ipack(domain%size(_FABM_VECTORIZED_DIMENSION_INDEX_)))
#  endif
#else
      n = 1
#endif

      call allocate_and_fill(cache%read, fill_values%read, 1, n)
      call allocate_and_fill(cache%read_hz, fill_values%read_hz, 1, n)
      call allocate_and_fill(cache%read_scalar, fill_values%read_scalar, 1, 1)
      call allocate_and_fill(cache%write, fill_values%write, 0, n)
   end subroutine

   subroutine create_horizontal_cache(domain, fill_values, cache)
      type (type_domain),            intent(in)  :: domain
      type (type_cache_fill_values), intent(in)  :: fill_values
      type (type_horizontal_cache),  intent(out) :: cache

      integer :: n, n_mod, i

#ifdef _HORIZONTAL_IS_VECTORIZED_
      n = domain%size(_FABM_VECTORIZED_DIMENSION_INDEX_)
      n_mod = mod(n, array_block_size)
      if (n_mod /= 0) n = n - n_mod + array_block_size
#  ifdef _HAS_MASK_
      allocate(cache%ipack(domain%size(_FABM_VECTORIZED_DIMENSION_INDEX_)))
#  endif
#else
      n = 1
#endif

      call allocate_and_fill(cache%read, fill_values%read, 1, n)
      call allocate_and_fill(cache%read_hz, fill_values%read_hz, 1, n)
      call allocate_and_fill(cache%read_scalar, fill_values%read_scalar, 1, 1)
      call allocate_and_fill(cache%write_hz, fill_values%write_hz, 0, n)
   end subroutine

   subroutine create_vertical_cache(domain, fill_values, cache)
      type (type_domain),            intent(in)  :: domain
      type (type_cache_fill_values), intent(in)  :: fill_values
      type (type_vertical_cache),    intent(out) :: cache

      integer :: n, n_mod, i

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      n = domain%size(_FABM_DEPTH_DIMENSION_INDEX_)
      n_mod = mod(n, array_block_size)
      if (n_mod /= 0) n = n - n_mod + array_block_size
#  ifdef _HAS_MASK_
      allocate(cache%ipack(domain%size(_FABM_DEPTH_DIMENSION_INDEX_)))
#  endif
#else
      n = 1
#endif

      call allocate_and_fill(cache%read, fill_values%read, 1, n)
      call allocate_and_fill(cache%read_hz, fill_values%read_hz, 1, 1)
      call allocate_and_fill(cache%read_scalar, fill_values%read_scalar, 1, 1)
      call allocate_and_fill(cache%write, fill_values%write, 0, n)
      call allocate_and_fill(cache%write_hz, fill_values%write_hz, 0, 1)
   end subroutine

subroutine begin_interior_task(domain, catalog, fill_values, task, cache _POSTARG_INTERIOR_IN_)
   type (type_domain),            intent(in)    :: domain
   type (type_catalog),           intent(in)    :: catalog
   type (type_cache_fill_values), intent(in)    :: fill_values
   type (type_task),              intent(in)    :: task
   type (type_interior_cache),    intent(inout) :: cache
   _DECLARE_ARGUMENTS_INTERIOR_IN_
   _DECLARE_INTERIOR_INDICES_

   integer :: i, j

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
#  ifdef _HAS_MASK_
   i = 0
   do _I_ = _START_, _STOP_
#    ifdef _FABM_HORIZONTAL_MASK_
      if (_IS_UNMASKED_(domain%mask_hz _INDEX_GLOBAL_HORIZONTAL_(_I_))) then
#    else
      if (_IS_UNMASKED_(domain%mask _INDEX_GLOBAL_INTERIOR_(_I_))) then
#    endif
          i = i + 1
          cache%ipack(i) = _I_
      end if
   end do
   _N_ = i
#  else
   _N_ = _STOP_ - _START_ + 1
#  endif
#endif
   
   do i = 1, size(task%load)
      j = task%load(i)
      if (j /= 0) then
         _PACK_GLOBAL_(catalog%interior(j)%p, cache%read, i, cache)
      end if
   end do
   do i = 1, size(task%load_hz)
      j = task%load_hz(i)
      if (j /= 0) then
         _HORIZONTAL_PACK_GLOBAL_(catalog%horizontal(j)%p, cache%read_hz, i, cache)
      end if
   end do
   do i = 1, size(task%load_scalar)
      j = task%load_scalar(i)
      if (j /= 0) then
         cache%read_scalar(i) = catalog%scalar(j)%p
      end if
   end do

   do i = 1, size(task%prefill)
      j = task%prefill(i)
      if (j == prefill_constant) then
         _CONCURRENT_LOOP_BEGIN_
            cache%write _INDEX_SLICE_PLUS_1_(i) = fill_values%write(i)
         _LOOP_END_
      elseif (j /= prefill_none) then
         _PACK_GLOBAL_(catalog%interior(j)%p, cache%write, i, cache)
      end if
   end do
end subroutine begin_interior_task

subroutine begin_horizontal_task(domain, catalog, fill_values, task, cache _POSTARG_HORIZONTAL_IN_)
   type (type_domain),            intent(in)    :: domain
   type (type_catalog),           intent(in)    :: catalog
   type (type_cache_fill_values), intent(in)    :: fill_values
   type (type_task),              intent(in)    :: task
   type (type_horizontal_cache),  intent(inout) :: cache
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   _DECLARE_HORIZONTAL_INDICES_

   integer :: i, j

#ifdef _HORIZONTAL_IS_VECTORIZED_
#  ifdef _HAS_MASK_
   i = 0
   do _J_ = _START_, _STOP_
      if (_IS_UNMASKED_(domain%mask_hz _INDEX_GLOBAL_HORIZONTAL_(_J_))) then
          i = i + 1
          cache%ipack(i) = _J_
      end if
   end do
   _N_ = i
#  else
   _N_ = _STOP_ - _START_ + 1
#  endif
#endif

   do i = 1, size(task%load_hz)
      j = task%load_hz(i)
      if (j /= 0) then
         _HORIZONTAL_PACK_GLOBAL_(catalog%horizontal(j)%p, cache%read_hz, i, cache)
      end if
   end do
   do i = 1, size(task%load_scalar)
      j = task%load_scalar(i)
      if (j /= 0) cache%read_scalar(i) = catalog%scalar(j)%p
   end do

   do i = 1, size(task%prefill_hz)
      if (task%prefill_hz(i) == prefill_constant) then
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) = fill_values%write_hz(i)
         _HORIZONTAL_LOOP_END_
      elseif (task%prefill_hz(i) /= prefill_none) then
         _HORIZONTAL_PACK_GLOBAL_(catalog%horizontal(task%prefill_hz(i))%p, cache%write_hz, i, cache)
      end if
   end do

   ! Also load boundary values for interior fields if performing surface or bottom-specific operations.
   if (task%operation == source_do_surface) then
      call load_surface_data(domain, catalog, task, cache _POSTARG_HORIZONTAL_IN_)
   elseif (task%operation == source_do_bottom) then
      call load_bottom_data (domain, catalog, task, cache _POSTARG_HORIZONTAL_IN_)
   end if
end subroutine begin_horizontal_task

subroutine load_surface_data(domain, catalog, task, cache _POSTARG_HORIZONTAL_IN_)
   type (type_domain),           intent(in)    :: domain
   type (type_catalog),          intent(in)    :: catalog
   type (type_task),             intent(in)    :: task
   type (type_horizontal_cache), intent(inout) :: cache
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   _DECLARE_HORIZONTAL_INDICES_

   integer :: i, j

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
   integer :: _VERTICAL_ITERATOR_
   _VERTICAL_ITERATOR_ = domain%surface_index
#endif

   do i = 1, size(task%load)
      j = task%load(i)
      if (j /= 0) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
#  ifdef _HAS_MASK_
            cache%read _INDEX_SLICE_PLUS_1_(i) = catalog%interior(j)%p _INDEX_GLOBAL_INTERIOR_(cache%ipack(_J_))
#  else
            cache%read _INDEX_SLICE_PLUS_1_(i) = catalog%interior(j)%p _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1)
#  endif
         _HORIZONTAL_LOOP_END_
#elif defined(_INTERIOR_IS_VECTORIZED_)
         cache%read(1,i) = catalog%interior(j)%p _INDEX_LOCATION_
#else
         cache%read(i) = catalog%interior(j)%p _INDEX_LOCATION_
#endif
      end if
   end do
end subroutine load_surface_data

subroutine load_bottom_data(domain, catalog, task, cache _POSTARG_HORIZONTAL_IN_)
   type (type_domain),           intent(in)    :: domain
   type (type_catalog),          intent(in)    :: catalog
   type (type_task),             intent(in)    :: task
   type (type_horizontal_cache), intent(inout) :: cache
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   _DECLARE_HORIZONTAL_INDICES_

   integer :: i, k
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
   integer :: _VERTICAL_ITERATOR_
#endif

#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_BOTTOM_INDEX_==0
   _VERTICAL_ITERATOR_ = domain%bottom_index
#endif

   do i = 1, size(task%load)
      k = task%load(i)
      if (k /= 0) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
#  if _FABM_BOTTOM_INDEX_==-1
#    ifdef _HAS_MASK_
            _VERTICAL_ITERATOR_ = domain%bottom_indices _INDEX_GLOBAL_HORIZONTAL_(cache%ipack(_J_))
#    else
            _VERTICAL_ITERATOR_ = domain%bottom_indices _INDEX_GLOBAL_HORIZONTAL_(_START_+_J_-1)
#    endif
#  endif
#  ifdef _HAS_MASK_
            cache%read _INDEX_SLICE_PLUS_1_(i) = catalog%interior(k)%p _INDEX_GLOBAL_INTERIOR_(cache%ipack(_J_))
#  else
            cache%read _INDEX_SLICE_PLUS_1_(i) = catalog%interior(k)%p _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1)
#  endif
         _HORIZONTAL_LOOP_END_
#elif defined(_INTERIOR_IS_VECTORIZED_)
#  if _FABM_BOTTOM_INDEX_==-1
         _VERTICAL_ITERATOR_ = domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#  endif
         cache%read(1,i) = catalog%interior(k)%p _INDEX_LOCATION_
#else
#  if _FABM_BOTTOM_INDEX_==-1
         _VERTICAL_ITERATOR_ = domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#  endif
         cache%read(i) = catalog%interior(k)%p _INDEX_LOCATION_
#endif
      end if
   end do
end subroutine load_bottom_data

subroutine begin_vertical_task(domain, catalog, fill_values, task, cache _POSTARG_VERTICAL_IN_)
   type (type_domain),            intent(in)    :: domain
   type (type_catalog),           intent(in)    :: catalog
   type (type_cache_fill_values), intent(in)    :: fill_values
   type (type_task),              intent(in)    :: task
   type (type_vertical_cache),    intent(inout) :: cache
   _DECLARE_ARGUMENTS_VERTICAL_IN_
   _DECLARE_VERTICAL_INDICES_

   integer :: i, j

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  ifdef _HAS_MASK_
   i = 0
   do _I_ = _VERTICAL_START_, _VERTICAL_STOP_
#    ifdef _FABM_HORIZONTAL_MASK_
      if (_IS_UNMASKED_(domain%mask_hz _INDEX_HORIZONTAL_LOCATION_)) then
#    else
      if ( _IS_UNMASKED_(domain%mask _INDEX_GLOBAL_VERTICAL_(_I_))) then
#    endif
          i = i + 1
          cache%ipack(i) = _I_
      end if
    end do
    _N_ = i
#  else
   _N_ = _VERTICAL_STOP_ - _VERTICAL_START_ + 1
#  endif
#endif

   do i = 1, size(task%load)
      j = task%load(i)
      if (j /= 0) then
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
         _CONCURRENT_VERTICAL_LOOP_BEGIN_
#  ifdef _HAS_MASK_
            cache%read _INDEX_SLICE_PLUS_1_(i) = catalog%interior(j)%p _INDEX_GLOBAL_VERTICAL_(cache%ipack(_I_))
#  else
            cache%read _INDEX_SLICE_PLUS_1_(i) = catalog%interior(j)%p _INDEX_GLOBAL_VERTICAL_(_VERTICAL_START_+_I_-1)
#  endif
         _VERTICAL_LOOP_END_
#elif defined(_INTERIOR_IS_VECTORIZED_)
         cache%read(1,i) = catalog%interior(j)%p _INDEX_LOCATION_
#else
         cache%read(i) = catalog%interior(j)%p _INDEX_LOCATION_
#endif
      end if
   end do

   do i = 1, size(task%load_hz)
      j = task%load_hz(i)
      if (j /= 0) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
         cache%read_hz(1, i) = catalog%horizontal(j)%p _INDEX_HORIZONTAL_LOCATION_
#else
         cache%read_hz(i) = catalog%horizontal(j)%p _INDEX_HORIZONTAL_LOCATION_
#endif
      end if
   end do

   do i = 1, size(task%load_scalar)
      j = task%load_scalar(i)
      if (j /= 0) cache%read_scalar(i) = catalog%scalar(j)%p
   end do

   do i = 1, size(task%prefill)
      if (task%prefill(i) == prefill_constant) then
#if defined(_INTERIOR_IS_VECTORIZED_)
         cache%write(:, i) = fill_values%write(i)
#else
         cache%write(i) = fill_values%write(i)
#endif
      elseif (task%prefill(i) /= prefill_none) then
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
         _CONCURRENT_VERTICAL_LOOP_BEGIN_
#  ifdef _HAS_MASK_
            cache%write _INDEX_SLICE_PLUS_1_(i) = catalog%interior(task%prefill(i))%p _INDEX_GLOBAL_VERTICAL_(cache%ipack(_I_))
#  else
            cache%write _INDEX_SLICE_PLUS_1_(i) = catalog%interior(task%prefill(i))%p _INDEX_GLOBAL_VERTICAL_(_VERTICAL_START_+_I_-1)
#  endif
         _VERTICAL_LOOP_END_
#elif defined(_INTERIOR_IS_VECTORIZED_)
         cache%write(1, i) = catalog%interior(task%prefill(i))%p _INDEX_LOCATION_
#else
         cache%write(i) = catalog%interior(task%prefill(i))%p _INDEX_LOCATION_
#endif
      end if
   end do

   do i = 1, size(task%prefill_hz)
      if (task%prefill_hz(i) == prefill_constant) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
         cache%write_hz(1, i) = fill_values%write_hz(i)
#else
         cache%write_hz(i) = fill_values%write_hz(i)
#endif
      elseif (task%prefill_hz(i) /= prefill_none) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
         cache%write_hz(1, i) = catalog%horizontal(task%prefill_hz(i))%p _INDEX_HORIZONTAL_LOCATION_
#else
         cache%write_hz(i) = catalog%horizontal(task%prefill_hz(i))%p _INDEX_HORIZONTAL_LOCATION_
#endif
      end if
   end do
end subroutine begin_vertical_task

subroutine end_interior_task(task, cache, store _POSTARG_INTERIOR_IN_)
   type (type_task),           intent(in)    :: task
   type (type_interior_cache), intent(in)    :: cache
   type (type_store),          intent(inout) :: store
   _DECLARE_ARGUMENTS_INTERIOR_IN_
   _DECLARE_INTERIOR_INDICES_

   integer :: i

   ! Copy newly written diagnostics that need to be saved to global store.
   do i = 1, size(task%save_sources)
      if (task%save_sources(i) /= 0) then
         _UNPACK_TO_GLOBAL_PLUS_1_(cache%write, task%save_sources(i), store%interior, i, cache, store%interior_missing_value(i))
      end if
   end do
end subroutine end_interior_task

subroutine end_horizontal_task(task, cache, store _POSTARG_HORIZONTAL_IN_)
   type (type_task),             intent(in)    :: task
   type (type_horizontal_cache), intent(in)    :: cache
   type (type_store),            intent(inout) :: store
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   _DECLARE_HORIZONTAL_INDICES_

   integer :: i

   ! Copy newly written horizontal diagnostics that need to be saved to global store.
   do i = 1, size(task%save_sources_hz)
      if (task%save_sources_hz(i) /= 0) then
         _HORIZONTAL_UNPACK_TO_GLOBAL_PLUS_1_(cache%write_hz, task%save_sources_hz(i), store%horizontal, i, cache, store%horizontal_missing_value(i))
      end if
   end do
end subroutine end_horizontal_task

subroutine end_vertical_task(task, cache, store _POSTARG_VERTICAL_IN_)
   type (type_task),           intent(in)    :: task
   type (type_vertical_cache), intent(in)    :: cache
   type (type_store),          intent(inout) :: store
   _DECLARE_ARGUMENTS_VERTICAL_IN_
   _DECLARE_VERTICAL_INDICES_

   integer :: i

   ! Copy diagnostics that need to be saved to global store.
   do i = 1, size(task%save_sources)
      if (task%save_sources(i) /= 0) then
         _VERTICAL_UNPACK_TO_GLOBAL_PLUS_1_(cache%write, task%save_sources(i), store%interior, i, cache, store%interior_missing_value(i))
      end if
   end do
   do i = 1, size(task%save_sources_hz)
      if (task%save_sources_hz(i) /= 0) then
         if (_N_ > 0) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
            store%horizontal(_PREARG_HORIZONTAL_LOCATION_ i) = cache%write_hz(1, task%save_sources_hz(i))
#else
            store%horizontal(_PREARG_HORIZONTAL_LOCATION_ i) = cache%write_hz(task%save_sources_hz(i))
#endif
         else
            store%horizontal(_PREARG_HORIZONTAL_LOCATION_ i) = store%horizontal_missing_value(i)
         end if
      end if
   end do
end subroutine end_vertical_task

   subroutine process_interior_slice(task, domain, catalog, cache_fill_values, store, cache _POSTARG_INTERIOR_IN_)
      type (type_task),              intent(in)    :: task
      type (type_domain),            intent(in)    :: domain
      type (type_catalog),           intent(in)    :: catalog
      type (type_cache_fill_values), intent(in)    :: cache_fill_values
      type (type_store),             intent(inout) :: store
      type (type_interior_cache),    intent(inout) :: cache
      _DECLARE_ARGUMENTS_INTERIOR_IN_

      type (type_call), pointer :: call_node
      integer                   :: i, j, k
      _DECLARE_INTERIOR_INDICES_

      call cache_pack(domain, catalog, cache_fill_values, task, cache _POSTARG_INTERIOR_IN_)

      call_node => task%first_call
      do while (associated(call_node))
         if (call_node%active) then
#ifndef NDEBUG
            call invalidate_interior_call_output(call_node, cache)
#endif

            select case (call_node%source)
            case (source_do);                    call call_node%model%do(cache)
            case (source_get_light_extinction);  call call_node%model%get_light_extinction(cache)
            case (source_get_vertical_movement); call call_node%model%get_vertical_movement(cache)
            end select

#ifndef NDEBUG
            call check_interior_call_output(call_node, cache)
#endif
         end if

         ! Copy outputs of interest to read cache so consecutive models can use it.
         _DO_CONCURRENT_(i,1,size(call_node%copy_commands_int))
            j = call_node%copy_commands_int(i)%read_index
            k = call_node%copy_commands_int(i)%write_index
            _CONCURRENT_LOOP_BEGIN_EX_(cache)
               cache%read _INDEX_SLICE_PLUS_1_(j) = cache%write _INDEX_SLICE_PLUS_1_(k)
            _LOOP_END_
         end do

         ! Move to next model
         call_node => call_node%next
      end do

      call cache_unpack(task, cache, store _POSTARG_INTERIOR_IN_)

   end subroutine process_interior_slice

   subroutine process_horizontal_slice(task, domain, catalog, cache_fill_values, store, cache _POSTARG_HORIZONTAL_IN_)
      type (type_task),              intent(in)    :: task
      type (type_domain),            intent(in)    :: domain
      type (type_catalog),           intent(in)    :: catalog
      type (type_cache_fill_values), intent(in)    :: cache_fill_values
      type (type_store),             intent(inout) :: store
      type (type_horizontal_cache),  intent(inout) :: cache
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_

      type (type_call), pointer :: call_node
      integer                   :: i, j, k
      _DECLARE_HORIZONTAL_INDICES_

      call cache_pack(domain, catalog, cache_fill_values, task, cache _POSTARG_HORIZONTAL_IN_)

      call_node => task%first_call
      do while (associated(call_node))
         if (call_node%active) then
#ifndef NDEBUG
            call invalidate_horizontal_call_output(call_node, cache)
#endif

            select case (call_node%source)
            case (source_do_surface);    call call_node%model%do_surface   (cache)
            case (source_do_bottom);     call call_node%model%do_bottom    (cache)
            case (source_do_horizontal); call call_node%model%do_horizontal(cache)
            end select

#ifndef NDEBUG
            call check_horizontal_call_output(call_node, cache)
#endif
         end if

         ! Copy outputs of interest to read cache so consecutive models can use it.
         _DO_CONCURRENT_(i,1,size(call_node%copy_commands_hz))
            j = call_node%copy_commands_hz(i)%read_index
            k = call_node%copy_commands_hz(i)%write_index
            _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache)
               cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
            _HORIZONTAL_LOOP_END_
         end do

         call_node => call_node%next
      end do

      call cache_unpack(task, cache, store _POSTARG_HORIZONTAL_IN_)

   end subroutine process_horizontal_slice

   subroutine process_vertical_slice(task, domain, catalog, cache_fill_values, store, cache _POSTARG_VERTICAL_IN_)
      type (type_task),              intent(in)    :: task
      type (type_domain),            intent(in)    :: domain
      type (type_catalog),           intent(in)    :: catalog
      type (type_cache_fill_values), intent(in)    :: cache_fill_values
      type (type_store),             intent(inout) :: store
      type (type_vertical_cache),    intent(inout) :: cache
      _DECLARE_ARGUMENTS_VERTICAL_IN_

      type (type_call), pointer :: call_node
      integer                   :: i, j, k
      _DECLARE_VERTICAL_INDICES_

      call cache_pack(domain, catalog, cache_fill_values, task, cache _POSTARG_VERTICAL_IN_)

      call_node => task%first_call
      do while (associated(call_node))

         if (call_node%active) then
#ifndef NDEBUG
            call invalidate_vertical_call_output(call_node, cache)
#endif

            call call_node%model%get_light(cache)

#ifndef NDEBUG
            call check_vertical_call_output(call_node, cache)
#endif
         end if

         ! Copy outputs of interest to read cache so consecutive models can use it.
         _DO_CONCURRENT_(i,1,size(call_node%copy_commands_int))
            j = call_node%copy_commands_int(i)%read_index
            k = call_node%copy_commands_int(i)%write_index
            _CONCURRENT_VERTICAL_LOOP_BEGIN_EX_(cache)
               cache%read _INDEX_SLICE_PLUS_1_(j) = cache%write _INDEX_SLICE_PLUS_1_(k)
            _VERTICAL_LOOP_END_
         end do
         _DO_CONCURRENT_(i,1,size(call_node%copy_commands_hz))
            j = call_node%copy_commands_hz(i)%read_index
            k = call_node%copy_commands_hz(i)%write_index
#ifdef _HORIZONTAL_IS_VECTORIZED_
            cache%read_hz(1, j) = cache%write_hz(1, k)
#else
            cache%read_hz(j) = cache%write_hz(k)
#endif
         end do

         call_node => call_node%next
      end do

      call cache_unpack(task, cache, store _POSTARG_VERTICAL_IN_)

   end subroutine process_vertical_slice

   subroutine invalidate_interior_call_output(call_node, cache)
      use fabm_graph, only: type_output_variable_set_node

      type (type_call),           intent(in)    :: call_node
      type (type_interior_cache), intent(inout) :: cache
      _DECLARE_INTERIOR_INDICES_

      type (type_output_variable_set_node), pointer :: output_variable
      integer                                       :: i

      output_variable => call_node%graph_node%outputs%first
      do while (associated(output_variable))
         if (output_variable%p%target%prefill == prefill_none) then
            i = output_variable%p%target%write_indices%value
#if defined(_INTERIOR_IS_VECTORIZED_)
            cache%write(:, i) = not_written
#else
            cache%write(i) = not_written
#endif
         end if
         output_variable => output_variable%next
      end do
   end subroutine

   subroutine invalidate_horizontal_call_output(call_node, cache)
      use fabm_graph, only: type_output_variable_set_node

      type (type_call),             intent(in)    :: call_node
      type (type_horizontal_cache), intent(inout) :: cache
      _DECLARE_INTERIOR_INDICES_

      type (type_output_variable_set_node), pointer :: output_variable
      integer                                       :: i

      output_variable => call_node%graph_node%outputs%first
      do while (associated(output_variable))
         if (output_variable%p%target%prefill == prefill_none) then
            i = output_variable%p%target%write_indices%value
#if defined(_HORIZONTAL_IS_VECTORIZED_)
            cache%write_hz(:, i) = not_written
#else
            cache%write_hz(i) = not_written
#endif
         end if
         output_variable => output_variable%next
      end do
   end subroutine

   subroutine invalidate_vertical_call_output(call_node, cache)
      use fabm_graph, only: type_output_variable_set_node

      type (type_call),           intent(in)    :: call_node
      type (type_vertical_cache), intent(inout) :: cache
      _DECLARE_INTERIOR_INDICES_

      type (type_output_variable_set_node), pointer :: output_variable
      integer                                       :: i

      output_variable => call_node%graph_node%outputs%first
      do while (associated(output_variable))
         if (output_variable%p%target%prefill == prefill_none) then
            i = output_variable%p%target%write_indices%value
            select case (output_variable%p%target%domain)
            case (domain_interior)
#if defined(_INTERIOR_IS_VECTORIZED_)
               cache%write(:, i) = not_written
#else
               cache%write(i) = not_written
#endif
            case (domain_surface, domain_bottom, domain_horizontal)
#if defined(_HORIZONTAL_IS_VECTORIZED_)
               cache%write_hz(:, i) = not_written
#else
               cache%write_hz(i) = not_written
#endif
            end select
         end if
         output_variable => output_variable%next
      end do
   end subroutine

   subroutine check_interior_call_output(call_node, cache)
      use fabm_graph, only: type_output_variable_set_node
      use, intrinsic :: ieee_arithmetic

      type (type_call),           intent(in) :: call_node
      type (type_interior_cache), intent(in) :: cache
      _DECLARE_INTERIOR_INDICES_

      type (type_output_variable_set_node), pointer :: output_variable

      output_variable => call_node%graph_node%outputs%first
      do while (associated(output_variable))
         _ASSERT_(output_variable%p%target%domain == domain_interior, 'check_interior_call_output', 'output not for interior domain')
         _LOOP_BEGIN_
            if (.not. ieee_is_finite(cache%write _INDEX_SLICE_PLUS_1_(output_variable%p%target%write_indices%value))) &
               call driver%fatal_error('check_interior_call_output', trim(call_node%model%get_path()) // ':' // trim(source2string(call_node%source)) // ' wrote non-finite data for ' // trim(output_variable%p%target%name))
         _LOOP_END_
         if (output_variable%p%target%prefill == prefill_none) then
            _LOOP_BEGIN_
               if (cache%write _INDEX_SLICE_PLUS_1_(output_variable%p%target%write_indices%value) == not_written) &
                  call driver%fatal_error('check_interior_call_output', trim(call_node%model%get_path()) // ':' // trim(source2string(call_node%source)) // ' failed to write data for ' // trim(output_variable%p%target%name))
            _LOOP_END_
         end if
         output_variable => output_variable%next
      end do
   end subroutine check_interior_call_output

   subroutine check_horizontal_call_output(call_node, cache)
      use fabm_graph, only: type_output_variable_set_node
      use, intrinsic :: ieee_arithmetic

      type (type_call),             intent(in) :: call_node
      type (type_horizontal_cache), intent(in) :: cache
      _DECLARE_HORIZONTAL_INDICES_

      type (type_output_variable_set_node), pointer :: output_variable

      output_variable => call_node%graph_node%outputs%first
      do while (associated(output_variable))
         _ASSERT_(iand(output_variable%p%target%domain, domain_horizontal) /= 0, 'check_horizontal_call_output', 'output not for horizontal domain')
         _HORIZONTAL_LOOP_BEGIN_
            if (.not. ieee_is_finite(cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(output_variable%p%target%write_indices%value))) &
               call driver%fatal_error('check_horizontal_call_output', trim(call_node%model%get_path()) // ':' // trim(source2string(call_node%source)) // ' wrote non-finite data for ' // trim(output_variable%p%target%name))
         _HORIZONTAL_LOOP_END_
         if (output_variable%p%target%prefill == prefill_none) then
            _HORIZONTAL_LOOP_BEGIN_
               if (cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(output_variable%p%target%write_indices%value) == not_written) &
                  call driver%fatal_error('check_horizontal_call_output', trim(call_node%model%get_path()) // ':' // trim(source2string(call_node%source)) // ' failed to write data for ' // trim(output_variable%p%target%name))
            _HORIZONTAL_LOOP_END_
         end if
         output_variable => output_variable%next
      end do
   end subroutine check_horizontal_call_output

   subroutine check_vertical_call_output(call_node, cache)
      use fabm_graph, only: type_output_variable_set_node
      use, intrinsic :: ieee_arithmetic

      type (type_call),           intent(in) :: call_node
      type (type_vertical_cache), intent(in) :: cache
      _DECLARE_VERTICAL_INDICES_

      type (type_output_variable_set_node), pointer :: output_variable

      output_variable => call_node%graph_node%outputs%first
      do while (associated(output_variable))
         select case (output_variable%p%target%domain)
         case (domain_interior)
            _VERTICAL_LOOP_BEGIN_
               if (.not. ieee_is_finite(cache%write _INDEX_SLICE_PLUS_1_(output_variable%p%target%write_indices%value))) &
                  call driver%fatal_error('check_vertical_call_output', trim(call_node%model%get_path()) // ':' // trim(source2string(call_node%source)) // ' wrote non-finite data for ' // trim(output_variable%p%target%name))
            _VERTICAL_LOOP_END_
            if (output_variable%p%target%prefill == prefill_none) then
               _VERTICAL_LOOP_BEGIN_
                  if (cache%write _INDEX_SLICE_PLUS_1_(output_variable%p%target%write_indices%value) == not_written) &
                     call driver%fatal_error('check_vertical_call_output', trim(call_node%model%get_path()) // ':' // trim(source2string(call_node%source)) // ' failed to write data for ' // trim(output_variable%p%target%name))
               _VERTICAL_LOOP_END_
            end if
         case (domain_surface, domain_bottom, domain_horizontal)
            if (.not. ieee_is_finite(cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(output_variable%p%target%write_indices%value))) &
               call driver%fatal_error('check_vertical_call_output', trim(call_node%model%get_path()) // ':' // trim(source2string(call_node%source)) // ' wrote non-finite data for ' // trim(output_variable%p%target%name))
            if (output_variable%p%target%prefill == prefill_none) then
               if (cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(output_variable%p%target%write_indices%value) == not_written) &
                  call driver%fatal_error('check_vertical_call_output', trim(call_node%model%get_path()) // ':' // trim(source2string(call_node%source)) // ' failed to write data for ' // trim(output_variable%p%target%name))
            end if
         end select
         output_variable => output_variable%next
      end do
   end subroutine check_vertical_call_output
   
end module
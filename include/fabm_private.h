#if _FABM_DIMENSION_COUNT_==0

! ---------------------
! 0D spatial context
! ---------------------

#define _LOCATION_
#define _LOCATION_DIMENSIONS_
#define _LOCATION_RANGE_

#elif _FABM_DIMENSION_COUNT_==1

! ---------------------
! 1D spatial context
! ---------------------

#define _LOCATION_ i__
#define _LOCATION_DIMENSIONS_ :
#define _LOCATION_RANGE_ istart__,istop__

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  define _GLOBAL_VERTICAL_(it) it
#else
#  define _BEGIN_OUTER_VERTICAL_LOOP_ do i__=istart__,istop__
#  define _END_OUTER_VERTICAL_LOOP_ end do
#endif

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
#  define _INTERIOR_FIXED_LOCATION_
#  define _GLOBAL_INTERIOR_(it) it
#else
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do i__=istart__,istop__
#  define _END_OUTER_INTERIOR_LOOP_ end do
#endif

#elif _FABM_DIMENSION_COUNT_==2

! ---------------------
! 2D spatial context
! ---------------------

#define _LOCATION_ i__,j__
#define _LOCATION_DIMENSIONS_ :,:
#define _LOCATION_RANGE_ istart__,istop__,jstart__,jstop__

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _HORIZONTAL_LOCATION_ j__
#    define _HORIZONTAL_LOCATION_RANGE_ jstart__,jstop__
#    define _BEGIN_OUTER_VERTICAL_LOOP_ do j__=jstart__,jstop__
#    define _END_OUTER_VERTICAL_LOOP_ end do
#    define _GLOBAL_VERTICAL_(it) it,j__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _HORIZONTAL_LOCATION_ i__
#    define _HORIZONTAL_LOCATION_RANGE_ istart__,istop__
#    define _BEGIN_OUTER_VERTICAL_LOOP_ do i__=istart__,istop__
#    define _END_OUTER_VERTICAL_LOOP_ end do
#    define _GLOBAL_VERTICAL_(it) i__,it
#  endif
#  define _HORIZONTAL_LOCATION_DIMENSIONS_ :
#else
#  define _BEGIN_OUTER_VERTICAL_LOOP_ do j__=jstart__,jstop__;do i__=istart__,istop__
#  define _END_OUTER_VERTICAL_LOOP_ end do;end do
#endif

#if _FABM_VECTORIZED_DIMENSION_INDEX_==1
#  define _INTERIOR_FIXED_LOCATION_ j__
#  define _GLOBAL_INTERIOR_(it) it,j__
#  if _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _GLOBAL_HORIZONTAL_(it) it
#  endif
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do j__=jstart__,jstop__
#  define _END_OUTER_INTERIOR_LOOP_ end do
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==2
#  define _INTERIOR_FIXED_LOCATION_ i__
#  define _GLOBAL_INTERIOR_(it) i__,it
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _GLOBAL_HORIZONTAL_(it) it
#  endif
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do i__=istart__,istop__
#  define _END_OUTER_INTERIOR_LOOP_ end do
#else
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do j__=jstart__,jstop__;do i__=istart__,istop__
#  define _END_OUTER_INTERIOR_LOOP_ end do;end do
#endif

#elif _FABM_DIMENSION_COUNT_==3

! ---------------------
! 3D spatial context
! ---------------------

#define _LOCATION_ i__,j__,k__
#define _LOCATION_DIMENSIONS_ :,:,:
#define _LOCATION_RANGE_ istart__,istop__,jstart__,jstop__,kstart__,kstop__

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _HORIZONTAL_LOCATION_ j__,k__
#    define _HORIZONTAL_LOCATION_RANGE_ jstart__,jstop__,kstart__,kstop__
#    define _BEGIN_OUTER_VERTICAL_LOOP_ do k__=kstart__,kstop__;do j__=jstart__,jstop__
#    define _END_OUTER_VERTICAL_LOOP_ end do;end do
#    define _GLOBAL_VERTICAL_(it) it,j__,k__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _HORIZONTAL_LOCATION_ i__,k__
#    define _HORIZONTAL_LOCATION_RANGE_ istart__,istop__,kstart__,kstop__
#    define _BEGIN_OUTER_VERTICAL_LOOP_ do k__=kstart__,kstop__;do i__=istart__,istop__
#    define _END_OUTER_VERTICAL_LOOP_ end do;end do
#    define _GLOBAL_VERTICAL_(it) i__,it,k__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _HORIZONTAL_LOCATION_ i__,j__
#    define _HORIZONTAL_LOCATION_RANGE_ istart__,istop__,jstart__,jstop__
#    define _BEGIN_OUTER_VERTICAL_LOOP_ do j__=jstart__,jstop__;do i__=istart__,istop__
#    define _END_OUTER_VERTICAL_LOOP_ end do;end do
#    define _GLOBAL_VERTICAL_(it) i__,j__,it
#  endif
#  define _HORIZONTAL_LOCATION_DIMENSIONS_ :,:
#else
#  define _BEGIN_OUTER_VERTICAL_LOOP_ do k__=kstart__,kstop__;do j__=jstart__,jstop__;do i__=istart__,istop__
#  define _END_OUTER_VERTICAL_LOOP_ end do;end do;end do
#endif

#if _FABM_VECTORIZED_DIMENSION_INDEX_==1
#  define _INTERIOR_FIXED_LOCATION_ j__,k__
#  define _GLOBAL_INTERIOR_(it) it,j__,k__
#  if _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _HORIZONTAL_FIXED_LOCATION_ k__
#    define _GLOBAL_HORIZONTAL_(it) it,k__
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do k__=kstart__,kstop__
#    define _END_OUTER_HORIZONTAL_LOOP_ end do
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _HORIZONTAL_FIXED_LOCATION_ j__
#    define _GLOBAL_HORIZONTAL_(it) it,j__
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do j__=jstart__,jstop__
#    define _END_OUTER_HORIZONTAL_LOOP_ end do
#  endif
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do k__=kstart__,kstop__;do j__=jstart__,jstop__
#  define _END_OUTER_INTERIOR_LOOP_ end do;end do
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==2
#  define _INTERIOR_FIXED_LOCATION_ i__,k__
#  define _GLOBAL_INTERIOR_(it) i__,it,k__
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _HORIZONTAL_FIXED_LOCATION_ k__
#    define _GLOBAL_HORIZONTAL_(it) it,k__
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do k__=kstart__,kstop__
#    define _END_OUTER_HORIZONTAL_LOOP_ end do
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _HORIZONTAL_FIXED_LOCATION_ i__
#    define _GLOBAL_HORIZONTAL_(it) i__,it
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do i__=istart__,istop__
#    define _END_OUTER_HORIZONTAL_LOOP_ end do
#  endif
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do k__=kstart__,kstop__;do i__=istart__,istop__
#  define _END_OUTER_INTERIOR_LOOP_ end do;end do
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==3
#  define _INTERIOR_FIXED_LOCATION_ i__,j__
#  define _GLOBAL_INTERIOR_(it) i__,j__,it
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _HORIZONTAL_FIXED_LOCATION_ j__
#    define _GLOBAL_HORIZONTAL_(it) j__,it
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do j__=jstart__,jstop__
#    define _END_OUTER_HORIZONTAL_LOOP_ end do
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _HORIZONTAL_FIXED_LOCATION_ i__
#    define _GLOBAL_HORIZONTAL_(it) i__,it
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do i__=istart__,istop__
#    define _END_OUTER_HORIZONTAL_LOOP_ end do
#  endif
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do j__=jstart__,jstop__;do i__=istart__,istop__
#  define _END_OUTER_INTERIOR_LOOP_ end do;end do
#else
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do k__=kstart__,kstop__;do j__=jstart__,jstop__;do i__=istart__,istop__
#  define _END_OUTER_INTERIOR_LOOP_ end do;end do;end do
#endif

#endif

#if _FABM_VECTORIZED_DIMENSION_INDEX_==1
#  define _ITERATOR_ i__
#  define _START_ istart__
#  define _STOP_ istop__
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==2
#  define _ITERATOR_ j__
#  define _START_ jstart__
#  define _STOP_ jstop__
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==3
#  define _ITERATOR_ k__
#  define _START_ kstart__
#  define _STOP_ kstop__
#endif

! If there is no depth dimension, horizontal dimensions match full dimensions.
#ifndef _FABM_DEPTH_DIMENSION_INDEX_
#  define _HORIZONTAL_FIXED_LOCATION_ _INTERIOR_FIXED_LOCATION_
#  define _HORIZONTAL_LOCATION_ _LOCATION_
#  define _HORIZONTAL_LOCATION_RANGE_ _LOCATION_RANGE_
#  define _HORIZONTAL_LOCATION_DIMENSIONS_ _LOCATION_DIMENSIONS_
#  define _HORIZONTAL_DIMENSION_COUNT_ _FABM_DIMENSION_COUNT_
#else
#  define _HORIZONTAL_DIMENSION_COUNT_ _FABM_DIMENSION_COUNT_-1
#endif

#if (!defined(_FABM_DEPTH_DIMENSION_INDEX_)||_FABM_DEPTH_DIMENSION_INDEX_==_FABM_VECTORIZED_DIMENSION_INDEX_)
#  define _BEGIN_OUTER_HORIZONTAL_LOOP_ _BEGIN_OUTER_INTERIOR_LOOP_
#  define _END_OUTER_HORIZONTAL_LOOP_ _END_OUTER_INTERIOR_LOOP_
#endif

#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&!defined(_FABM_DEPTH_DIMENSION_INDEX_)
#  define _GLOBAL_HORIZONTAL_(it) _GLOBAL_INTERIOR_(it)
#endif

! Check for additional required preprocessor variables.
#ifndef _LOCATION_
#  error BUG: Preprocessor variable _LOCATION_ must be defined.
#endif
#ifndef _LOCATION_DIMENSIONS_
#  error BUG: Preprocessor variable _LOCATION_DIMENSIONS_ must be defined.
#endif
#if _HORIZONTAL_DIMENSION_COUNT_>0
#  ifndef _HORIZONTAL_LOCATION_
#    error BUG: Preprocessor variable _HORIZONTAL_LOCATION_ must be defined.
#  endif
#  ifndef _LOCATION_DIMENSIONS_
#    error BUG: Preprocessor variable _HORIZONTAL_LOCATION_DIMENSIONS_ must be defined.
#  endif
#endif
#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&!defined(_GLOBAL_INTERIOR_)
#  error BUG: Preprocessor variable _GLOBAL_INTERIOR_ must be defined since _FABM_VECTORIZED_DIMENSION_INDEX_ is set.
#endif
#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&!defined(_GLOBAL_VERTICAL_)
#  error BUG: Preprocessor variable _GLOBAL_VERTICAL_ must be defined since _FABM_DEPTH_DIMENSION_INDEX_ is set.
#endif
#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&_FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_&&!defined(_GLOBAL_HORIZONTAL_)
#  error BUG: Preprocessor variable _GLOBAL_HORIZONTAL_ must be defined since _FABM_VECTORIZED_DIMENSION_INDEX_ is set and not equal to _FABM_DEPTH_DIMENSION_INDEX_.
#endif

! =======================================================================================================
! Process spatial mask, based on the following variables provided by the driver:
!   _FABM_MASK_TYPE_ (data type of mask elements, e.g., logical, integer or real)
!   _FABM_MASKED_VALUE_ or _FABM_UNMASKED_VALUE_ (mask value for masked and unmasked cells, respectively)
! =======================================================================================================

#ifdef _FABM_MASK_TYPE_
#  define _HAS_MASK_
#endif

#ifdef _HAS_MASK_
#  ifndef _INTERIOR_IS_VECTORIZED_
#    error _FABM_MASK_TYPE_/_FABM_MASKED_VALUE_/_FABM_UNMASKED_VALUE_ are not used if no dimension is vectorized.
#  endif
#  ifdef _FABM_IS_UNMASKED_
#    define _IS_UNMASKED_(maskvalue) _FABM_IS_UNMASKED_(maskvalue)
#  elif defined(_FABM_MASKED_VALUE_)
#    define _IS_UNMASKED_(maskvalue) maskvalue/=_FABM_MASKED_VALUE_
#  elif defined(_FABM_UNMASKED_VALUE_)
#    define _IS_UNMASKED_(maskvalue) maskvalue==_FABM_UNMASKED_VALUE_
#  else
#    error If _FABM_MASK_TYPE_ is set, _FABM_MASKED_VALUE_, _FABM_UNMASKED_VALUE_ or _FABM_IS_UNMASKED_ must be set as well.
#  endif
#else
#  ifdef _FABM_IS_UNMASKED_
#    error To use _FABM_IS_UNMASKED_, _FABM_MASK_TYPE_ must be set as well.
#  endif
#  ifdef _FABM_MASKED_VALUE_
#    error To use _FABM_MASKED_VALUE_, _FABM_MASK_TYPE_ must be set as well.
#  endif
#  ifdef _FABM_UNMASKED_VALUE_
#    error To use _FABM_UNMASKED_VALUE_, _FABM_MASK_TYPE_ must be set as well.
#  endif
#  define _IS_UNMASKED_(maskvalue) .true.
#endif

#ifdef _FABM_CONTIGUOUS_
#  define _CONTIGUOUS_ ,contiguous
#else
#  define _CONTIGUOUS_
#endif

#if _FABM_DIMENSION_COUNT_==0||(_FABM_DIMENSION_COUNT_==1&&defined(_FABM_VECTORIZED_DIMENSION_INDEX_))
#  define _BEGIN_OUTER_INTERIOR_LOOP_
#  define _END_OUTER_INTERIOR_LOOP_
#endif

#ifndef _BEGIN_OUTER_HORIZONTAL_LOOP_
#  define _BEGIN_OUTER_HORIZONTAL_LOOP_
#  define _END_OUTER_HORIZONTAL_LOOP_
#endif

#if _FABM_DIMENSION_COUNT_==0||(_FABM_DIMENSION_COUNT_==1&&defined(_FABM_DEPTH_DIMENSION_INDEX_))
#  define _BEGIN_OUTER_VERTICAL_LOOP_
#  define _END_OUTER_VERTICAL_LOOP_
#endif

! =================================================================================
! Further preprocessor macros for specifying spatial dimensionality and position
! =================================================================================


#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
#  define _DIMENSION_EXT_SLICE_ ,dimension(_START_:)
#  define _DIMENSION_EXT_SLICE_PLUS_1_ ,dimension(_START_:,:)
#  define _DIMENSION_EXT_SLICE_PLUS_2_ ,dimension(_START_:,:,:)
#  define _INDEX_EXT_SLICE_ (_START_+_I_-1)
#  define _INDEX_EXT_SLICE_PLUS_1_(i) (_START_+_I_-1,i)
#  define _INDEX_EXT_SLICE_PLUS_2_(i,j) (_START_+_I_-1,i,j)
#else
#  define _DIMENSION_EXT_SLICE_
#  define _DIMENSION_EXT_SLICE_PLUS_1_ ,dimension(:)
#  define _DIMENSION_EXT_SLICE_PLUS_2_ ,dimension(:,:)
#  define _INDEX_EXT_SLICE_
#  define _INDEX_EXT_SLICE_PLUS_1_(i) (i)
#  define _INDEX_EXT_SLICE_PLUS_2_(i,j) (i,j)
#endif

#ifdef _HORIZONTAL_IS_VECTORIZED_
! Horizontal fields are 1D
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_ _DIMENSION_EXT_SLICE_
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_ _DIMENSION_EXT_SLICE_PLUS_1_
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_2_ _DIMENSION_EXT_SLICE_PLUS_2_
#  define _INDEX_EXT_HORIZONTAL_SLICE_ _INDEX_EXT_SLICE_
#  define _INDEX_EXT_HORIZONTAL_SLICE_PLUS_1_(i) _INDEX_EXT_SLICE_PLUS_1_(i)
#  define _INDEX_EXT_HORIZONTAL_SLICE_PLUS_2_(i,j) _INDEX_EXT_SLICE_PLUS_2_(i,j)
#else
! Horizontal fields are 0D
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_ ,dimension(:)
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_2_ ,dimension(:,:)
#  define _INDEX_EXT_HORIZONTAL_SLICE_
#  define _INDEX_EXT_HORIZONTAL_SLICE_PLUS_1_(i) (i)
#  define _INDEX_EXT_HORIZONTAL_SLICE_PLUS_2_(i,j) (i,j)
#endif

! ---------------------------------------------------------------------------------
! Dimension attribute and index specifyer for horizontal (2D) fields.
! ---------------------------------------------------------------------------------

#if _HORIZONTAL_DIMENSION_COUNT_>0
#  define _INDEX_HORIZONTAL_LOCATION_ (_HORIZONTAL_LOCATION_)
#  define _DIMENSION_GLOBAL_HORIZONTAL_ ,dimension(_HORIZONTAL_LOCATION_DIMENSIONS_)
#  define _ATTRIBUTES_GLOBAL_HORIZONTAL_ _DIMENSION_GLOBAL_HORIZONTAL_ _CONTIGUOUS_
#  define _ARG_HORIZONTAL_LOCATION_ _HORIZONTAL_LOCATION_
#  define _ARGUMENTS_HORIZONTAL_LOCATION_ ,_ARG_HORIZONTAL_LOCATION_
#  define _ARGUMENTS_HORIZONTAL_LOCATION_RANGE_ ,_HORIZONTAL_LOCATION_RANGE_
#  define _PREARG_HORIZONTAL_LOCATION_ _ARG_HORIZONTAL_LOCATION_,
#  define _PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ _HORIZONTAL_LOCATION_DIMENSIONS_,
#  define _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_ integer,intent(in) :: _HORIZONTAL_LOCATION_
#  define _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_RANGE_ integer,intent(in) :: _HORIZONTAL_LOCATION_RANGE_
#else
#  define _INDEX_HORIZONTAL_LOCATION_
#  define _DIMENSION_GLOBAL_HORIZONTAL_
#  define _ATTRIBUTES_GLOBAL_HORIZONTAL_
#  define _ARG_HORIZONTAL_LOCATION_
#  define _ARGUMENTS_HORIZONTAL_LOCATION_
#  define _ARGUMENTS_HORIZONTAL_LOCATION_RANGE_
#  define _PREARG_HORIZONTAL_LOCATION_
#  define _PREARG_HORIZONTAL_LOCATION_DIMENSIONS_
#  define _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_
#  define _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_RANGE_
#endif

! ---------------------------------------------------------------------------------
! Dimension attribute and index specifyer for full 3D fields.
! ---------------------------------------------------------------------------------

#if _FABM_DIMENSION_COUNT_>0
#  define _INDEX_LOCATION_ (_LOCATION_)
#  define _DIMENSION_GLOBAL_ ,dimension(_LOCATION_DIMENSIONS_)
#  define _ATTRIBUTES_GLOBAL_ _DIMENSION_GLOBAL_ _CONTIGUOUS_
#  define _POSTARG_LOCATION_ ,_LOCATION_
#  define _ARGUMENTS_LOCATION_RANGE_ ,_LOCATION_RANGE_
#  define _PREARG_LOCATION_ _LOCATION_,
#  define _PREARG_LOCATION_DIMENSIONS_ _LOCATION_DIMENSIONS_,
#  define _DECLARE_ARGUMENTS_LOCATION_ integer,intent(in) :: _LOCATION_
#  define _DECLARE_ARGUMENTS_LOCATION_RANGE_ integer,intent(in) :: _LOCATION_RANGE_
#  define _DECLARE_LOCATION_ integer :: _LOCATION_
#else
#  define _INDEX_LOCATION_
#  define _DIMENSION_GLOBAL_
#  define _ATTRIBUTES_GLOBAL_
#  define _POSTARG_LOCATION_
#  define _ARGUMENTS_LOCATION_RANGE_
#  define _PREARG_LOCATION_
#  define _PREARG_LOCATION_DIMENSIONS_
#  define _DECLARE_ARGUMENTS_LOCATION_
#  define _DECLARE_ARGUMENTS_LOCATION_RANGE_
#  define _DECLARE_LOCATION_
#endif

#define _DIMENSION_GLOBAL_PLUS_1_ ,dimension(_PREARG_LOCATION_DIMENSIONS_ :)
#define _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ ,dimension(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ :)

#ifdef _GLOBAL_INTERIOR_
!  Interior is vectorized; forward provided iterator
#  define _INDEX_GLOBAL_INTERIOR_(it) (_GLOBAL_INTERIOR_(it))
#  define _INDEX_GLOBAL_INTERIOR_PLUS_1_(it,j) (_GLOBAL_INTERIOR_(it),j)
#else
!  Interior is not vectorized; just index to local point in space.
#  define _INDEX_GLOBAL_INTERIOR_(it) _INDEX_LOCATION_
#  define _INDEX_GLOBAL_INTERIOR_PLUS_1_(it,j) (_PREARG_LOCATION_ j)
#endif

#ifdef _GLOBAL_HORIZONTAL_
!  Interior is vectorized; forward provided iterator
#  define _INDEX_GLOBAL_HORIZONTAL_(it) (_GLOBAL_HORIZONTAL_(it))
#  define _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(it,j) (_GLOBAL_HORIZONTAL_(it),j)
#else
!  Interior is not vectorized; just index to local point in space.
#  define _INDEX_GLOBAL_HORIZONTAL_(it) _INDEX_HORIZONTAL_LOCATION_
#  define _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(it,j) (_PREARG_HORIZONTAL_LOCATION_ j)
#endif

#ifdef _GLOBAL_VERTICAL_
!  Interior is vectorized; forward provided iterator
#  define _INDEX_GLOBAL_VERTICAL_(it) (_GLOBAL_VERTICAL_(it))
#  define _INDEX_GLOBAL_VERTICAL_PLUS_1_(it,j) (_GLOBAL_VERTICAL_(it),j)
#else
!  Interior is not vectorized; just index to local point in space.
#  define _INDEX_GLOBAL_VERTICAL_(it) _INDEX_LOCATION_
#  define _INDEX_GLOBAL_VERTICAL_PLUS_1_(it,j) (_PREARG_LOCATION_ j)
#endif

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
!  ---------------------------------------------------------------------------------
!  INTERIOR procedures operate on a data slice over one spatial dimension.
!  ---------------------------------------------------------------------------------
#  if _FABM_DIMENSION_COUNT_>1
#    define _ARG_INTERIOR_FIXED_LOCATION_ ,_INTERIOR_FIXED_LOCATION_
#  else
#    define _ARG_INTERIOR_FIXED_LOCATION_
#  endif
#  define _ARG_INTERIOR_IN_ _START_,_STOP_ _ARG_INTERIOR_FIXED_LOCATION_
#  define _POSTARG_INTERIOR_IN_ ,_ARG_INTERIOR_IN_
#  define _PREARG_INTERIOR_IN_ _ARG_INTERIOR_IN_,
#  define _DECLARE_ARGUMENTS_INTERIOR_IN_ integer,intent(in) :: _START_,_STOP_ _ARG_INTERIOR_FIXED_LOCATION_
#else
!  ---------------------------------------------------------------------------------
!  INTERIOR procedures operate on one point at a time.
!  ---------------------------------------------------------------------------------
#  define _ARG_INTERIOR_IN_ _LOCATION_
#  define _POSTARG_INTERIOR_IN_ _POSTARG_LOCATION_
#  define _PREARG_INTERIOR_IN_ _PREARG_LOCATION_
#  define _DECLARE_ARGUMENTS_INTERIOR_IN_ _DECLARE_ARGUMENTS_LOCATION_
#endif

#ifdef _HORIZONTAL_IS_VECTORIZED_
!  ---------------------------------------------------------------------------------
!  HORIZONTAL procedures operate on a data slice over one spatial dimension.
!  This will be the same dimension that INTERIOR procedures operate upon.
!  ---------------------------------------------------------------------------------
#  if (_FABM_DIMENSION_COUNT_>2||(_FABM_DIMENSION_COUNT_==2&&!defined(_FABM_DEPTH_DIMENSION_INDEX_)))
#    define _ARG_HORIZONTAL_FIXED_LOCATION_ ,_HORIZONTAL_FIXED_LOCATION_
#  else
#    define _ARG_HORIZONTAL_FIXED_LOCATION_
#  endif
#  define _ARG_HORIZONTAL_IN_ _START_,_STOP_ _ARG_HORIZONTAL_FIXED_LOCATION_
#  define _POSTARG_HORIZONTAL_IN_ ,_ARG_HORIZONTAL_IN_
#  define _PREARG_HORIZONTAL_IN_ _ARG_HORIZONTAL_IN_,
#  define _DECLARE_ARGUMENTS_HORIZONTAL_IN_ integer,intent(in) :: _START_,_STOP_ _ARG_HORIZONTAL_FIXED_LOCATION_
#else
!  ---------------------------------------------------------------------------------
!  HORIZONTAL procedures operate on one point at a time.
!  ---------------------------------------------------------------------------------
#  define _ARG_HORIZONTAL_IN_ _ARG_HORIZONTAL_LOCATION_
#  define _POSTARG_HORIZONTAL_IN_ _ARGUMENTS_HORIZONTAL_LOCATION_
#  define _PREARG_HORIZONTAL_IN_ _PREARG_HORIZONTAL_LOCATION_
#  define _DECLARE_ARGUMENTS_HORIZONTAL_IN_ _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_
#endif

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
!  ---------------------------------------------------------------------------------
!  VERTICAL procedures operate on a data slice over one spatial dimension.
!  ---------------------------------------------------------------------------------
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _VERTICAL_ITERATOR_ i__
#    define _VERTICAL_START_ istart__
#    define _VERTICAL_STOP_ istop__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _VERTICAL_ITERATOR_ j__
#    define _VERTICAL_START_ jstart__
#    define _VERTICAL_STOP_ jstop__
#  else
#    define _VERTICAL_ITERATOR_ k__
#    define _VERTICAL_START_ kstart__
#    define _VERTICAL_STOP_ kstop__
#  endif
#  if _FABM_DIMENSION_COUNT_==1
#    define _ARG_VERTICAL_FIXED_LOCATION_
#  else
#    define _ARG_VERTICAL_FIXED_LOCATION_ ,_HORIZONTAL_LOCATION_
#  endif
#  define _ARG_VERTICAL_IN_ _VERTICAL_START_,_VERTICAL_STOP_ _ARG_VERTICAL_FIXED_LOCATION_
#  define _POSTARG_VERTICAL_IN_ ,_ARG_VERTICAL_IN_
#  define _PREARG_VERTICAL_IN_ _ARG_VERTICAL_IN_,
#  define _DECLARE_ARGUMENTS_VERTICAL_IN_ integer,intent(in) :: _VERTICAL_START_,_VERTICAL_STOP_ _ARG_VERTICAL_FIXED_LOCATION_
#else
!  ---------------------------------------------------------------------------------
!  VERTICAL procedures operate on one point at a time.
!  ---------------------------------------------------------------------------------
#  define _ARG_VERTICAL_IN_ _LOCATION_
#  define _POSTARG_VERTICAL_IN_ _POSTARG_LOCATION_
#  define _PREARG_VERTICAL_IN_ _PREARG_LOCATION_
#  define _DECLARE_ARGUMENTS_VERTICAL_IN_ _DECLARE_ARGUMENTS_LOCATION_
#endif

#ifdef _HAS_MASK_
#  define _PACK_GLOBAL_(in,out,i,cache) _CONCURRENT_LOOP_BEGIN_EX_(cache);out _INDEX_SLICE_PLUS_1_(i) = in _INDEX_GLOBAL_INTERIOR_(cache%ipack(_I_));_LOOP_END_
#  define _PACK_GLOBAL_PLUS_1_(in,i,out,j,cache) _CONCURRENT_LOOP_BEGIN_EX_(cache);out _INDEX_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_INTERIOR_PLUS_1_(cache%ipack(_I_),i);_LOOP_END_
#  define _UNPACK_(in,i,out,cache,missing) _DO_CONCURRENT_(_I_,_START_,_STOP_);out(_I_) = in(cache%iunpack(_I_),i);end do
#  define _UNPACK_TO_PLUS_1_(in,i,out,j,cache,missing) _DO_CONCURRENT_(_I_,_START_,_STOP_);out(_I_,j) = in(cache%iunpack(_I_),i);end do
#  define _UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j,cache) _DO_CONCURRENT_(_I_,_START_,_STOP_);out(_I_,j) = out(_I_,j) + in(cache%iunpack(_I_),i);end do
#  define _UNPACK_TO_GLOBAL_(in,i,out,cache,missing) _DO_CONCURRENT_(_I_,_START_,_STOP_);out _INDEX_GLOBAL_INTERIOR_(_I_) = in(cache%iunpack(_I_),i);end do
#  define _UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,cache,missing) _DO_CONCURRENT_(_I_,_START_,_STOP_);out _INDEX_GLOBAL_INTERIOR_PLUS_1_(_I_,j) = in(cache%iunpack(_I_),i);end do
#else
#  define _PACK_GLOBAL_(in,out,i,cache) _CONCURRENT_LOOP_BEGIN_EX_(cache);out _INDEX_SLICE_PLUS_1_(i) = in _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1);_LOOP_END_
#  define _PACK_GLOBAL_PLUS_1_(in,i,out,j,cache) _CONCURRENT_LOOP_BEGIN_EX_(cache);out _INDEX_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_INTERIOR_PLUS_1_(_START_+_I_-1,i);_LOOP_END_
#  define _UNPACK_(in,i,out,cache,missing) _CONCURRENT_LOOP_BEGIN_EX_(cache);out _INDEX_EXT_SLICE_ = in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#  define _UNPACK_TO_PLUS_1_(in,i,out,j,cache,missing) _CONCURRENT_LOOP_BEGIN_EX_(cache);out _INDEX_EXT_SLICE_PLUS_1_(j) = in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#  define _UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j,cache) _CONCURRENT_LOOP_BEGIN_EX_(cache);out _INDEX_EXT_SLICE_PLUS_1_(j) = out _INDEX_EXT_SLICE_PLUS_1_(j) + in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#  define _UNPACK_TO_GLOBAL_(in,i,out,cache,missing) _CONCURRENT_LOOP_BEGIN_EX_(cache);out _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1) = in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#  define _UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,cache,missing) _CONCURRENT_LOOP_BEGIN_EX_(cache);out _INDEX_GLOBAL_INTERIOR_PLUS_1_(_START_+_I_-1,j) = in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#endif

#if defined(_HORIZONTAL_IS_VECTORIZED_)&&defined(_HAS_MASK_)
#  define _HORIZONTAL_PACK_GLOBAL_(in,out,j,cache) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache);out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_HORIZONTAL_(cache%ipack(_J_));_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_PACK_GLOBAL_PLUS_1_(in,i,out,j,cache) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache);out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(cache%ipack(_J_),i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_(in,i,out,cache,missing) _DO_CONCURRENT_(_J_,_START_,_STOP_);out(_J_) = in(cache%iunpack(_J_),i);end do
#  define _HORIZONTAL_UNPACK_TO_PLUS_1_(in,i,out,j,cache,missing) _DO_CONCURRENT_(_J_,_START_,_STOP_);out(_J_,j) = in(cache%iunpack(_J_),i);end do
#  define _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j,cache) _DO_CONCURRENT_(_J_,_START_,_STOP_);out(_J_,j) = out(_J_,j) + in(cache%iunpack(_J_),i);end do
#  define _HORIZONTAL_UNPACK_TO_GLOBAL_(in,i,out,cache,missing) _DO_CONCURRENT_(_J_,_START_,_STOP_);out _INDEX_GLOBAL_HORIZONTAL_(_J_) = in(cache%iunpack(_J_),i);end do
#  define _HORIZONTAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,cache,missing) _DO_CONCURRENT_(_J_,_START_,_STOP_);out _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_J_,j) = in(cache%iunpack(_J_),i);end do
#else
#  define _HORIZONTAL_PACK_GLOBAL_(in,out,j,cache) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache);out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_HORIZONTAL_(_START_+_J_-1);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_PACK_GLOBAL_PLUS_1_(in,i,out,j,cache) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache);out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_START_+_J_-1,i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_(in,i,out,cache,missing) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache);out _INDEX_EXT_HORIZONTAL_SLICE_ = in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_TO_PLUS_1_(in,i,out,j,cache,missing) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache);out _INDEX_EXT_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j,cache) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache);out _INDEX_EXT_HORIZONTAL_SLICE_PLUS_1_(j) = out _INDEX_EXT_HORIZONTAL_SLICE_PLUS_1_(j) + in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_TO_GLOBAL_(in,i,out,cache,missing) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache);out _INDEX_GLOBAL_HORIZONTAL_(_START_+_J_-1) = in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,cache,missing) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(cache);out _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_START_+_J_-1,j) = in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_HORIZONTAL_LOOP_END_
#endif

#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&defined(_HAS_MASK_)
#  define _VERTICAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,cache,missing) _DO_CONCURRENT_(_I_,_VERTICAL_START_,_VERTICAL_STOP_);out _INDEX_GLOBAL_VERTICAL_PLUS_1_(_I_,j) = in(cache%iunpack(_I_),i);end do
#else
#  define _VERTICAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,cache,missing) _CONCURRENT_VERTICAL_LOOP_BEGIN_EX_(cache);out _INDEX_GLOBAL_VERTICAL_PLUS_1_(_VERTICAL_START_+_I_-1,j) = in _INDEX_SLICE_PLUS_1_(i);_VERTICAL_LOOP_END_
#endif

#if _FABM_DIMENSION_COUNT_==0

! ---------------------
! 0D spatial context
! ---------------------

#define _LOCATION_
#define _LOCATION_DIMENSIONS_

#elif _FABM_DIMENSION_COUNT_==1

! ---------------------
! 1D spatial context
! ---------------------

#define _LOCATION_ i__
#define _LOCATION_DIMENSIONS_ :

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  define _GLOBAL_VERTICAL_(it) it
#endif

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
#  define _INTERIOR_FIXED_LOCATION_
#  define _GLOBAL_INTERIOR_(it) it
#endif

#elif _FABM_DIMENSION_COUNT_==2

! ---------------------
! 2D spatial context
! ---------------------

#define _LOCATION_ i__,j__
#define _LOCATION_DIMENSIONS_ :,:

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _HORIZONTAL_LOCATION_ j__
#    define _GLOBAL_VERTICAL_(it) it,j__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _HORIZONTAL_LOCATION_ i__
#    define _GLOBAL_VERTICAL_(it) i__,it
#  endif
#  define _HORIZONTAL_LOCATION_DIMENSIONS_ :
#endif

#if _FABM_VECTORIZED_DIMENSION_INDEX_==1
#  define _INTERIOR_FIXED_LOCATION_ j__
#  define _GLOBAL_INTERIOR_(it) it,j__
#  if _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _GLOBAL_HORIZONTAL_(it) it
#  endif
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==2
#  define _INTERIOR_FIXED_LOCATION_ i__
#  define _GLOBAL_INTERIOR_(it) i__,it
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _GLOBAL_HORIZONTAL_(it) it
#  endif
#endif

#elif _FABM_DIMENSION_COUNT_==3

! ---------------------
! 3D spatial context
! ---------------------

#define _LOCATION_ i__,j__,k__
#define _LOCATION_DIMENSIONS_ :,:,:

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _HORIZONTAL_LOCATION_ j__,k__
#    define _GLOBAL_VERTICAL_(it) it,j__,k__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _HORIZONTAL_LOCATION_ i__,k__
#    define _GLOBAL_VERTICAL_(it) i__,it,k__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _HORIZONTAL_LOCATION_ i__,j__
#    define _GLOBAL_VERTICAL_(it) i__,j__,it
#  endif
#  define _HORIZONTAL_LOCATION_DIMENSIONS_ :,:
#endif

#if _FABM_VECTORIZED_DIMENSION_INDEX_==1
#  define _INTERIOR_FIXED_LOCATION_ j__,k__
#  define _GLOBAL_INTERIOR_(it) it,j__,k__
#  if _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _HORIZONTAL_FIXED_LOCATION_ k__
#    define _GLOBAL_HORIZONTAL_(it) it,k__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _HORIZONTAL_FIXED_LOCATION_ j__
#    define _GLOBAL_HORIZONTAL_(it) it,j__
#  endif
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==2
#  define _INTERIOR_FIXED_LOCATION_ i__,k__
#  define _GLOBAL_INTERIOR_(it) i__,it,k__
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _HORIZONTAL_FIXED_LOCATION_ k__
#    define _GLOBAL_HORIZONTAL_(it) it,k__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _HORIZONTAL_FIXED_LOCATION_ i__
#    define _GLOBAL_HORIZONTAL_(it) i__,it
#  endif
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==3
#  define _INTERIOR_FIXED_LOCATION_ i__,j__
#  define _GLOBAL_INTERIOR_(it) i__,j__,it
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _HORIZONTAL_FIXED_LOCATION_ j__
#    define _GLOBAL_HORIZONTAL_(it) j__,it
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _HORIZONTAL_FIXED_LOCATION_ i__
#    define _GLOBAL_HORIZONTAL_(it) i__,it
#  endif
#endif

#endif

#if _FABM_DEPTH_DIMENSION_INDEX_==1
#  define _VERTICAL_ITERATOR_ i__
#elif _FABM_DEPTH_DIMENSION_INDEX_==2
#  define _VERTICAL_ITERATOR_ j__
#elif _FABM_DEPTH_DIMENSION_INDEX_==3
#  define _VERTICAL_ITERATOR_ k__
#endif

! If there is no depth dimension, horizontal dimensions match full dimensions.
#ifndef _FABM_DEPTH_DIMENSION_INDEX_
#  define _HORIZONTAL_FIXED_LOCATION_ _INTERIOR_FIXED_LOCATION_
#  define _HORIZONTAL_LOCATION_ _LOCATION_
#  define _HORIZONTAL_LOCATION_DIMENSIONS_ _LOCATION_DIMENSIONS_
#  define _HORIZONTAL_DIMENSION_COUNT_ _FABM_DIMENSION_COUNT_
#else
#  define _HORIZONTAL_DIMENSION_COUNT_ _FABM_DIMENSION_COUNT_-1
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
#endif

! =================================================================================
! Further preprocessor macros for specifying spatial dimensionality and position
! =================================================================================


#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
#  define _DIMENSION_EXT_SLICE_ ,dimension(loop_start:)
#  define _DIMENSION_EXT_SLICE_PLUS_1_ ,dimension(loop_start:,:)
#  define _DIMENSION_EXT_SLICE_PLUS_2_ ,dimension(loop_start:,:,:)
#  define _INDEX_EXT_SLICE_ (loop_start+_I_-1)
#  define _INDEX_EXT_SLICE_PLUS_1_(i) (loop_start+_I_-1,i)
#  define _INDEX_EXT_SLICE_PLUS_2_(i,j) (loop_start+_I_-1,i,j)
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
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_ ,dimension(loop_start:)
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_ ,dimension(loop_start:,:)
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_2_ ,dimension(loop_start:,:,:)
#else
! Horizontal fields are 0D
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_ ,dimension(:)
#  define _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_2_ ,dimension(:,:)
#endif

! ---------------------------------------------------------------------------------
! Dimension attribute and index specifyer for horizontal (2D) fields.
! ---------------------------------------------------------------------------------

#if _HORIZONTAL_DIMENSION_COUNT_>0
#  define _INDEX_HORIZONTAL_LOCATION_ (_HORIZONTAL_LOCATION_)
#  define _DIMENSION_GLOBAL_HORIZONTAL_ ,dimension(_HORIZONTAL_LOCATION_DIMENSIONS_)
#  define _ARGUMENTS_HORIZONTAL_LOCATION_ ,_HORIZONTAL_LOCATION_
#  define _PREARG_HORIZONTAL_LOCATION_ _HORIZONTAL_LOCATION_,
#  define _PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ _HORIZONTAL_LOCATION_DIMENSIONS_,
#  define _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_ integer,intent(in) :: _HORIZONTAL_LOCATION_
#else
#  define _INDEX_HORIZONTAL_LOCATION_
#  define _DIMENSION_GLOBAL_HORIZONTAL_
#  define _ARGUMENTS_HORIZONTAL_LOCATION_
#  define _PREARG_HORIZONTAL_LOCATION_
#  define _PREARG_HORIZONTAL_LOCATION_DIMENSIONS_
#  define _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_
#endif

! ---------------------------------------------------------------------------------
! Dimension attribute and index specifyer for full 3D fields.
! ---------------------------------------------------------------------------------

#if _FABM_DIMENSION_COUNT_>0
#  define _INDEX_LOCATION_ (_LOCATION_)
#  define _DIMENSION_GLOBAL_ ,dimension(_LOCATION_DIMENSIONS_)
#  define _ARGUMENTS_LOCATION_ ,_LOCATION_
#  define _PREARG_LOCATION_ _LOCATION_,
#  define _PREARG_LOCATION_DIMENSIONS_ _LOCATION_DIMENSIONS_,
#  define _DECLARE_ARGUMENTS_LOCATION_ integer,intent(in) :: _LOCATION_
#else
#  define _INDEX_LOCATION_
#  define _DIMENSION_GLOBAL_
#  define _ARGUMENTS_LOCATION_
#  define _PREARG_LOCATION_
#  define _PREARG_LOCATION_DIMENSIONS_
#  define _DECLARE_ARGUMENTS_LOCATION_
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
#  define _ARGUMENTS_INTERIOR_LENGTH_ ,_N_
#  define _ARGUMENTS_INTERIOR_IN_ ,loop_start,loop_stop _ARG_INTERIOR_FIXED_LOCATION_
#  define _DECLARE_ARGUMENTS_INTERIOR_IN_ integer,intent(in) :: loop_start,loop_stop _ARG_INTERIOR_FIXED_LOCATION_
#else
!  ---------------------------------------------------------------------------------
!  INTERIOR procedures operate on one point at a time.
!  ---------------------------------------------------------------------------------
#  define _ARGUMENTS_INTERIOR_LENGTH_
#  define _ARGUMENTS_INTERIOR_IN_ _ARGUMENTS_LOCATION_
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
#  define _ARGUMENTS_HORIZONTAL_LENGTH_ _ARGUMENTS_INTERIOR_LENGTH_
#  define _ARGUMENTS_HORIZONTAL_IN_ ,loop_start,loop_stop _ARG_HORIZONTAL_FIXED_LOCATION_
#  define _DECLARE_ARGUMENTS_HORIZONTAL_IN_ integer,intent(in) :: loop_start,loop_stop _ARG_HORIZONTAL_FIXED_LOCATION_

#else
!  ---------------------------------------------------------------------------------
!  HORIZONTAL procedures operate on one point at a time.
!  ---------------------------------------------------------------------------------
#  define _ARGUMENTS_HORIZONTAL_LENGTH_
#  define _ARGUMENTS_HORIZONTAL_IN_ _ARGUMENTS_HORIZONTAL_LOCATION_
#  define _DECLARE_ARGUMENTS_HORIZONTAL_IN_ _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_
#endif

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
!  ---------------------------------------------------------------------------------
!  VERTICAL procedures operate on a data slice over one spatial dimension.
!  ---------------------------------------------------------------------------------
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _VERTICAL_ITERATOR_ i__
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _VERTICAL_ITERATOR_ j__
#  else
#    define _VERTICAL_ITERATOR_ k__
#  endif
#  if _FABM_DIMENSION_COUNT_==1
#    define _ARG_VERTICAL_FIXED_LOCATION_
#  else
#    define _ARG_VERTICAL_FIXED_LOCATION_ ,_HORIZONTAL_LOCATION_
#  endif
#  define _ARGUMENTS_VERTICAL_LENGTH_ ,_N_
#  define _ARGUMENTS_VERTICAL_IN_ ,loop_start,loop_stop _ARG_VERTICAL_FIXED_LOCATION_
#  define _DECLARE_ARGUMENTS_VERTICAL_IN_ integer,intent(in) :: loop_start,loop_stop _ARG_VERTICAL_FIXED_LOCATION_
#else
!  ---------------------------------------------------------------------------------
!  VERTICAL procedures operate on one point at a time.
!  ---------------------------------------------------------------------------------
#  define _ARGUMENTS_VERTICAL_LENGTH_
#  define _ARGUMENTS_VERTICAL_IN_ _ARGUMENTS_LOCATION_
#  define _DECLARE_ARGUMENTS_VERTICAL_IN_ _DECLARE_ARGUMENTS_LOCATION_
#endif

#ifdef _HAS_MASK_
! Using pack/unpack intrinsics
!#  define _PACK_GLOBAL_(in, out, i, env) out(:, i) = pack(in _INDEX_GLOBAL_INTERIOR_(loop_start:loop_stop), env%mask)
!#  define _PACK_GLOBAL_PLUS_1_(in, i, out, j, env) out(:, j) = pack(in _INDEX_GLOBAL_INTERIOR_PLUS_1_(loop_start:loop_stop, i), env%mask)
!#  define _UNPACK_(in, i, out, env, missing) out(:) = unpack(in(:, i), env%mask, missing)
!#  define _UNPACK_TO_PLUS_1_(in, i, out, j, env, missing) out(:, j) = unpack(in(:, i), env%mask, missing)
!#  define _UNPACK_AND_ADD_TO_PLUS_1_(in, i, out, j, env) out(:, j) = out(:, j) + unpack(in(:, i), env%mask, 0._rk)
!#  define _UNPACK_TO_GLOBAL_(in, i, out, env, missing) out _INDEX_GLOBAL_INTERIOR_(loop_start:loop_stop) = unpack(in(:, i), env%mask, missing)
!#  define _UNPACK_TO_GLOBAL_PLUS_1_(in, i, out, j, env, missing) out _INDEX_GLOBAL_INTERIOR_PLUS_1_(loop_start:loop_stop, j) = unpack(in(:, i), env%mask, missing)

! Using our own arrays with pack/env indices
#  define _PACK_GLOBAL_(in,out,i,env) _CONCURRENT_LOOP_BEGIN_EX_(env);out _INDEX_SLICE_PLUS_1_(i) = in _INDEX_GLOBAL_INTERIOR_(env%ipack(_I_));_LOOP_END_
#  define _PACK_GLOBAL_PLUS_1_(in,i,out,j,env) _CONCURRENT_LOOP_BEGIN_;out _INDEX_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_INTERIOR_PLUS_1_(env%ipack(_I_),i);_LOOP_END_
#  define _UNPACK_(in,i,out,env,missing) _DO_CONCURRENT_(_I_,loop_start,loop_stop);if (env%iunpack(_I_)/=0) then;out(_I_) = in(env%iunpack(_I_),i);else;out(_I_) = missing;end if;_LOOP_END_
#  define _UNPACK_TO_PLUS_1_(in,i,out,j,env,missing) _DO_CONCURRENT_(_I_,loop_start,loop_stop);if (env%iunpack(_I_)/=0) then;out(_I_,j) = in(env%iunpack(_I_),i);else;out(_I_,j) = missing;end if;_LOOP_END_
!#  define _UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j,env) _DO_CONCURRENT_(_I_,loop_start,loop_stop);if (env%iunpack(_I_)/=0) then;out(_I_,j) = out(_I_,j) + in(env%iunpack(_I_),i);end if;_LOOP_END_
#  define _UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j,env) _CONCURRENT_LOOP_BEGIN_EX_(env);out(env%ipack(_I_),j) = out(env%ipack(_I_),j) + in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#  define _UNPACK_TO_GLOBAL_(in,i,out,env,missing) _DO_CONCURRENT_(_I_,loop_start,loop_stop);if (env%iunpack(_I_)/=0) then;out _INDEX_GLOBAL_INTERIOR_(_I_) = in(env%iunpack(_I_),i);else;out _INDEX_GLOBAL_INTERIOR_(_I_) = missing;end if;_LOOP_END_
#  define _UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,env,missing) _DO_CONCURRENT_(_I_,loop_start,loop_stop);if (env%iunpack(_I_)/=0) then;out _INDEX_GLOBAL_INTERIOR_PLUS_1_(_I_,j) = in(env%iunpack(_I_),i);else;out _INDEX_GLOBAL_INTERIOR_PLUS_1_(_I_,j) = missing;end if;_LOOP_END_
#else
#  define _PACK_GLOBAL_(in,out,i,env) _CONCURRENT_LOOP_BEGIN_EX_(env);out _INDEX_SLICE_PLUS_1_(i) = in _INDEX_GLOBAL_INTERIOR_(loop_start+_I_-1);_LOOP_END_
#  define _PACK_GLOBAL_PLUS_1_(in,i,out,j,env) _CONCURRENT_LOOP_BEGIN_EX_(env);out _INDEX_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_INTERIOR_PLUS_1_(loop_start+_I_-1,i);_LOOP_END_
#  define _UNPACK_(in,i,out,env,missing) _CONCURRENT_LOOP_BEGIN_EX_(env);out _INDEX_EXT_SLICE_ = in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#  define _UNPACK_TO_PLUS_1_(in,i,out,j,env,missing) _CONCURRENT_LOOP_BEGIN_EX_(env);out _INDEX_EXT_SLICE_PLUS_1_(j) = in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#  define _UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j,env) _CONCURRENT_LOOP_BEGIN_EX_(env);out _INDEX_EXT_SLICE_PLUS_1_(j) = out _INDEX_EXT_SLICE_PLUS_1_(j) + in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#  define _UNPACK_TO_GLOBAL_(in,i,out,env,missing) _CONCURRENT_LOOP_BEGIN_EX_(env);out _INDEX_GLOBAL_INTERIOR_(loop_start+_I_-1) = in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#  define _UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,env,missing) _CONCURRENT_LOOP_BEGIN_EX_(env);out _INDEX_GLOBAL_INTERIOR_PLUS_1_(loop_start+_I_-1,j) = in _INDEX_SLICE_PLUS_1_(i);_LOOP_END_
#endif

#if defined(_HORIZONTAL_IS_VECTORIZED_)&&defined(_HAS_MASK_)
! Using pack/unpack intrinsics
!#  define _HORIZONTAL_PACK_GLOBAL_(in,out,j,env) out(:,j) = pack(in _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop),env%mask)
!#  define _HORIZONTAL_PACK_GLOBAL_PLUS_1_(in,i,out,j, env) out(:,j) = pack(in _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(loop_start:loop_stop,i), env%mask)
!#  define _HORIZONTAL_UNPACK_TO_PLUS_1_(in,i,out,j, env,missing) out(:,j) = unpack(in(:,i), env%mask,missing)
!#  define _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j, env) out(:,j) = out(:,j) + unpack(in(:,i), env%mask,0._rk)
!#  define _HORIZONTAL_UNPACK_TO_GLOBAL_(in,i,out,env,missing) out _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop) = unpack(in(:,i), env%mask,missing)
!#  define _HORIZONTAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,env,missing) out _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(loop_start:loop_stop,j) = unpack(in(:,i),env%mask,missing)

! Using our own arrays with pack/unpack indices
#  define _HORIZONTAL_PACK_GLOBAL_(in,out,j,env) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(env);out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_HORIZONTAL_(env%ipack(_J_));_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_PACK_GLOBAL_PLUS_1_(in,i,out,j,env) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_;out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(env%ipack(_J_),i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_TO_PLUS_1_(in,i,out,j,env,missing) _DO_CONCURRENT_(_J_,loop_start,loop_stop);if (env%iunpack(_J_)/=0) then;out(_J_,j) = in(env%iunpack(_J_),i);else;out(_J_,j) = missing;end if;_LOOP_END_
!#  define _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j, env) _DO_CONCURRENT_(_J_,loop_start,loop_stop);if (env%iunpack(_J_)/=0) then;out(_J_,j) = out(_J_,j) + in(env%iunpack(_J_),i);end if;_LOOP_END_
#  define _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j,env) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(env);out(env%ipack(_J_),j) = out(env%ipack(_J_),j) + in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_LOOP_END_
#  define _HORIZONTAL_UNPACK_TO_GLOBAL_(in,i,out,env,missing) _DO_CONCURRENT_(_J_,loop_start,loop_stop);if (env%iunpack(_J_)/=0) then;out _INDEX_GLOBAL_HORIZONTAL_(_J_) = in(env%iunpack(_J_),i);else;out _INDEX_GLOBAL_HORIZONTAL_(_J_) = missing;end if;_LOOP_END_
#  define _HORIZONTAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,env,missing) _DO_CONCURRENT_(_J_,loop_start,loop_stop);if (env%iunpack(_J_)/=0) then;out _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_J_,j) = in(env%iunpack(_J_),i);else;out _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_J_,j) = missing;end if;_LOOP_END_
#else
#  define _HORIZONTAL_PACK_GLOBAL_(in,out,j,env) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(env);out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_HORIZONTAL_(loop_start+_J_-1);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_PACK_GLOBAL_PLUS_1_(in,i,out,j,env) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(env);out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(loop_start+_J_-1,i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_TO_PLUS_1_(in,i,out,j,env,missing) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(env);out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(in,i,out,j,env) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(env);out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = out _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) + in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_TO_GLOBAL_(in,i,out,env,missing) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(env);out _INDEX_GLOBAL_HORIZONTAL_(loop_start+_J_-1) = in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_HORIZONTAL_LOOP_END_
#  define _HORIZONTAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,env,missing) _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(env);out _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(loop_start+_J_-1,j) = in _INDEX_HORIZONTAL_SLICE_PLUS_1_(i);_HORIZONTAL_LOOP_END_
#endif

#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&defined(_HAS_MASK_)
! Using pack/unpack intrinsics
!#  define _VERTICAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,env,missing) out _INDEX_GLOBAL_VERTICAL_PLUS_1_(loop_start:loop_stop,j) = unpack(in(:,i),env%mask,missing)

! Using our own arrays with pack/unpack indices
#  define _VERTICAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,env,missing) _DO_CONCURRENT_(_I_,loop_start,loop_stop);if (env%iunpack(_I_)/=0) then;out _INDEX_GLOBAL_VERTICAL_PLUS_1_(_I_,j) = in(env%iunpack(_I_),i);else;out _INDEX_GLOBAL_VERTICAL_PLUS_1_(_I_,j) = missing;end if;_LOOP_END_
#else
#  define _VERTICAL_UNPACK_TO_GLOBAL_PLUS_1_(in,i,out,j,env,missing) _CONCURRENT_VERTICAL_LOOP_BEGIN_EX_(env);out _INDEX_GLOBAL_VERTICAL_PLUS_1_(loop_start+_I_-1,j) = in _INDEX_SLICE_PLUS_1_(i);_VERTICAL_LOOP_END_
#endif

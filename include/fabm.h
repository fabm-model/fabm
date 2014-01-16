! ========================================================
! Process spatial domain and set derived macros:
!   _LOCATION_
!	_LOCATION_DIMENSIONS_
!   _LOCATION_HZ_
!   _LOCATION_HZ_DIMENSIONS_
!   _VARIABLE_1DLOOP_
!   _LOCATION_1DLOOP_
!   _FABM_HORIZONTAL_IS_SCALAR_
! ========================================================

#ifndef _FABM_DIMENSION_COUNT_
#error Preprocessor variable _FABM_DIMENSION_COUNT_ must be defined.
#endif

#if (_FABM_DIMENSION_COUNT_<0||_FABM_DIMENSION_COUNT_>3)
#error Preprocessor variable _FABM_DIMENSION_COUNT_ takes values between 0 and 3 only.
#endif

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#if (_FABM_DEPTH_DIMENSION_INDEX_<1)||(_FABM_DEPTH_DIMENSION_INDEX_>_FABM_DIMENSION_COUNT_)
#error Preprocessor variable _FABM_DEPTH_DIMENSION_INDEX_ takes values between 1 and _FABM_DIMENSION_COUNT_ only.
#endif
#endif

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
#if (_FABM_VECTORIZED_DIMENSION_INDEX_<1)||(_FABM_VECTORIZED_DIMENSION_INDEX_>_FABM_DIMENSION_COUNT_)
#error Preprocessor variable _FABM_VECTORIZED_DIMENSION_INDEX_ takes values between 1 and _FABM_DIMENSION_COUNT_ only.
#endif
#endif

#if _FABM_DIMENSION_COUNT_==0

! ---------------------
! 0D spatial context
! ---------------------

#define _LOCATION_
#define _LOCATION_DIMENSIONS_
#define _FABM_HORIZONTAL_IS_SCALAR_

#elif _FABM_DIMENSION_COUNT_==1

! ---------------------
! 1D spatial context
! ---------------------

#define _LOCATION_ i__
#define _LOCATION_DIMENSIONS_ :

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#define _FABM_HORIZONTAL_IS_SCALAR_
#endif

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
#define _VARIABLE_1DLOOP_ i__
#define _LOCATION_1DLOOP_
#endif

#elif _FABM_DIMENSION_COUNT_==2

! ---------------------
! 2D spatial context
! ---------------------

#define _LOCATION_ i__,j__
#define _LOCATION_DIMENSIONS_ :,:

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#if _FABM_DEPTH_DIMENSION_INDEX_==1
#define _LOCATION_HZ_ j__
#elif _FABM_DEPTH_DIMENSION_INDEX_==2
#define _LOCATION_HZ_ i__
#endif
#define _LOCATION_DIMENSIONS_HZ_ :
#endif

#if _FABM_VECTORIZED_DIMENSION_INDEX_==1
#define _VARIABLE_1DLOOP_ i__
#define _LOCATION_1DLOOP_ j__
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==2
#define _VARIABLE_1DLOOP_ j__
#define _LOCATION_1DLOOP_ i__
#endif

#elif _FABM_DIMENSION_COUNT_==3

! ---------------------
! 3D spatial context
! ---------------------

#define _LOCATION_ i__,j__,k__
#define _LOCATION_DIMENSIONS_ :,:,:

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#if _FABM_DEPTH_DIMENSION_INDEX_==1
#define _LOCATION_HZ_ j__,k__
#elif _FABM_DEPTH_DIMENSION_INDEX_==2
#define _LOCATION_HZ_ i__,k__
#elif _FABM_DEPTH_DIMENSION_INDEX_==3
#define _LOCATION_HZ_ i__,j__
#endif
#define _LOCATION_DIMENSIONS_HZ_ :,:
#endif

#if _FABM_VECTORIZED_DIMENSION_INDEX_==1
#define _VARIABLE_1DLOOP_ i__
#define _LOCATION_1DLOOP_ j__,k__
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==2
#define _VARIABLE_1DLOOP_ j__
#define _LOCATION_1DLOOP_ i__,k__
#elif _FABM_VECTORIZED_DIMENSION_INDEX_==3
#define _VARIABLE_1DLOOP_ k__
#define _LOCATION_1DLOOP_ i__,j__
#endif

#endif

! If there is no depth dimension, horizontal dimensions match full dimensions.
#ifndef _FABM_DEPTH_DIMENSION_INDEX_
#define _LOCATION_HZ_ _LOCATION_
#define _LOCATION_DIMENSIONS_HZ_ _LOCATION_DIMENSIONS_
#endif

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
#define _FABM_USE_1D_LOOP_
#endif

#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&(_FABM_DEPTH_DIMENSION_INDEX_!=_FABM_VECTORIZED_DIMENSION_INDEX_)
#define _FABM_USE_1D_LOOP_IN_HORIZONTAL_
#endif

! Check for additional required preprocessor variables.
#ifndef _LOCATION_
#error Preprocessor variable _LOCATION_ must be defined.
#endif
#ifndef _LOCATION_DIMENSIONS_
#error Preprocessor variable _LOCATION_DIMENSIONS_ must be defined.
#endif
#ifndef _FABM_HORIZONTAL_IS_SCALAR_
#ifndef _LOCATION_HZ_
#error Preprocessor variable _LOCATION_HZ_ must be defined.
#endif
#ifndef _LOCATION_DIMENSIONS_
#error Preprocessor variable _LOCATION_DIMENSIONS_HZ_ must be defined.
#endif
#endif

! =======================================================================================================
! Process spatial mask, based on the following variables provided by the driver:
!   _FABM_MASK_TYPE_ (data type of mask elements, e.g., logical, integer or real)
!   _FABM_MASKED_VALUE_ or _FABM_UNMASKED_VALUE_ (mask value for masked and unmasked cells, respectively)
! =======================================================================================================

#ifdef _FABM_MASK_TYPE_
#define _FABM_MASK_
#endif

#ifdef _FABM_MASK_

#ifndef _FABM_USE_1D_LOOP_
#error _FABM_MASK_TYPE_/_FABM_MASKED_VALUE_/_FABM_UNMASKED_VALUE_ are not used if no dimension is vectorized.
#endif

#ifdef _FABM_MASKED_VALUE_
#define _FABM_IS_UNMASKED_(maskvalue) maskvalue/=_FABM_MASKED_VALUE_
#elif defined(_FABM_UNMASKED_VALUE_)
#define _FABM_IS_UNMASKED_(maskvalue) maskvalue==_FABM_UNMASKED_VALUE_
#else
#error If _FABM_MASK_TYPE_ is set, _FABM_MASKED_VALUE_ and/or _FABM_UNMASKED_VALUE_ must be set as well.
#endif

#endif

! =====================
! Portability macros
! =====================

#ifndef _FABM_REAL_KIND_
#define _FABM_REAL_KIND_ selected_real_kind(13)
#endif 

! Constants related to floating point precision; used throughout FABM.
#undef REALTYPE
#undef _ZERO_
#undef _ONE_
#define REALTYPE real(rk)
#define _ZERO_ 0._rk
#define _ONE_  1._rk

! For backward compatibility only [pre Fortran 2003]:
#define _CLASS_ class
#define _ALLOCATABLE_ allocatable
#define _NULL_
#define _ALLOCATED_ allocated

! =================================================================================
! Further preprocessor macros for specifying spatial dimensionality and position
! =================================================================================

! ---------------------------------------------------------------------------------
! Dimension attribute and index specifyer for horizontal (2D) fields.
! ---------------------------------------------------------------------------------

#ifdef _FABM_HORIZONTAL_IS_SCALAR_
#define _INDEX_LOCATION_HZ_
#define _DIMENSION_GLOBAL_HORIZONTAL_
#define _PREARG_LOCATION_HZ_
#define _PREARG_LOCATION_DIMENSIONS_HZ_
#else
#define _INDEX_LOCATION_HZ_ (_LOCATION_HZ_)
#define _DIMENSION_GLOBAL_HORIZONTAL_ ,dimension(_LOCATION_DIMENSIONS_HZ_)
#define _PREARG_LOCATION_HZ_ _LOCATION_HZ_,
#define _PREARG_LOCATION_DIMENSIONS_HZ_ _LOCATION_DIMENSIONS_HZ_,
#endif

! ---------------------------------------------------------------------------------
! Dimension attribute and index specifyer for full 3D fields.
! ---------------------------------------------------------------------------------

#if _FABM_DIMENSION_COUNT_>0
#define _INDEX_LOCATION_ (_LOCATION_)
#define _DIMENSION_GLOBAL_ ,dimension(_LOCATION_DIMENSIONS_)
#define _ARG_LOCATION_ ,_LOCATION_
#define _ARG_LOCATION_DIMENSIONS_ ,_LOCATION_DIMENSIONS_
#define _PREARG_LOCATION_ _LOCATION_,
#define _PREARG_LOCATION_DIMENSIONS_ _LOCATION_DIMENSIONS_,
#define _DECLARE_LOCATION_ARG_ integer,intent(in) :: _LOCATION_
#else
#define _INDEX_LOCATION_
#define _DIMENSION_GLOBAL_
#define _ARG_LOCATION_
#define _ARG_LOCATION_DIMENSIONS_
#define _PREARG_LOCATION_
#define _PREARG_LOCATION_DIMENSIONS_
#define _DECLARE_LOCATION_ARG_
#endif

#define _DIMENSION_GLOBAL_PLUS_1_ ,dimension(_PREARG_LOCATION_DIMENSIONS_ :)
#define _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ ,dimension(_PREARG_LOCATION_DIMENSIONS_HZ_ :)

#ifdef _FABM_USE_1D_LOOP_

! ---------------------------------------------------------------------------------
! 1D vectorized: FABM subroutines operate on one spatial dimension.
! ---------------------------------------------------------------------------------

#if _FABM_DIMENSION_COUNT_>1
#define _ARG_LOCATION_1DLOOP_ ,_LOCATION_1DLOOP_
#else
#define _ARG_LOCATION_1DLOOP_
#endif

! Length of the vectorized dimension, used in fabm_initialize only.
#define _DOMAIN_SIZE_1D_ _VARIABLE_1DLOOP_

! Dummy argument and argument declaration for location specification.
#define _ARG_LOCATION_ND_ ,fabm_loop_start,fabm_loop_stop _ARG_LOCATION_1DLOOP_
#define _DECLARE_LOCATION_ARG_ND_ integer,intent(in) :: fabm_loop_start,fabm_loop_stop _ARG_LOCATION_1DLOOP_;integer :: _VARIABLE_1DLOOP_

! Beginning and end of spatial loop
#ifdef _FABM_MASK_
#define _LOOP_BEGIN_EX_(environment) do _VARIABLE_1DLOOP_=fabm_loop_start,fabm_loop_stop;if (_FABM_IS_UNMASKED_(environment%mask _INDEX_LOCATION_)) then
#define _LOOP_END_ end if;end do
#else
#define _LOOP_BEGIN_EX_(environment) do _VARIABLE_1DLOOP_=fabm_loop_start,fabm_loop_stop
#define _LOOP_END_ end do
#endif

! Dimensionality of generic space-dependent arguments.
#define _DIMENSION_SLICE_ ,dimension(:)
#define _DIMENSION_SLICE_PLUS_1_ ,dimension(:,:)
#define _DIMENSION_SLICE_PLUS_2_ ,dimension(:,:,:)

#define _INDEX_OUTPUT_ (_VARIABLE_1DLOOP_-fabm_loop_start+1)
#define _INDEX_OUTPUT_1D_(index) (_VARIABLE_1DLOOP_-fabm_loop_start+1,index)
#define _INDEX_OUTPUT_2D_(index1,index2) (_VARIABLE_1DLOOP_-fabm_loop_start+1,index1,index2)

#define _SIZE_SLICE_ fabm_loop_stop-fabm_loop_start+1,
#define _DIMENSION_SLICE_AUTOMATIC_ ,dimension(fabm_loop_stop-fabm_loop_start+1)

#else

! ---------------------------------------------------------------------------------
! Not vectorized: FABM subroutines operate one the local state only.
! ---------------------------------------------------------------------------------

! Length of the vectorized dimension, used in fabm_initialize only.
#define _DOMAIN_SIZE_1D_ 1

! Dummy argument and argument declaration for location specification.
#define _ARG_LOCATION_ND_ _ARG_LOCATION_
#define _DECLARE_LOCATION_ARG_ND_ _DECLARE_LOCATION_ARG_

! Beginning and end of spatial loop
#define _LOOP_BEGIN_EX_(environment)
#define _LOOP_END_

! Dimensionality of generic space-dependent arguments.
#define _DIMENSION_SLICE_
#define _DIMENSION_SLICE_PLUS_1_ ,dimension(:)
#define _DIMENSION_SLICE_PLUS_2_ ,dimension(:,:)

! Expressions for indexing space-dependent FABM variables defined on the full spatial domain.
#define _INDEX_OUTPUT_
#define _INDEX_OUTPUT_1D_(index) (index)
#define _INDEX_OUTPUT_2D_(index1,index2) (index1,index2)

#define _SIZE_SLICE_
#define _DIMENSION_SLICE_AUTOMATIC_

#endif

#define _LOOP_BEGIN_ _LOOP_BEGIN_EX_(environment)

#ifdef _FABM_USE_1D_LOOP_IN_HORIZONTAL_

! ---------------------------------------------------------------------------------
! Functions operating on horizontonal slices of the domain will be vectorized just like the full domain.
! ---------------------------------------------------------------------------------

#define _ARG_LOCATION_VARS_HZ_ _ARG_LOCATION_ND_
#define _DECLARE_LOCATION_ARG_HZ_ _DECLARE_LOCATION_ARG_ND_

! Spatial loop for quantities defined on horizontal slice of the full spatial domain.
#define _HORIZONTAL_LOOP_BEGIN_EX_(environment) _LOOP_BEGIN_EX_(environment)
#define _HORIZONTAL_LOOP_END_ _LOOP_END_

! Vertical dimension is not among those vectorized:
! dimensionality of horizontal arrays will be equal to that of full domain arrays.
#define _DIMENSION_SLICE_HORIZONTAL_ _DIMENSION_SLICE_
#define _DIMENSION_SLICE_HORIZONTAL_PLUS_1_ _DIMENSION_SLICE_PLUS_1_
#define _DIMENSION_SLICE_HORIZONTAL_PLUS_2_ _DIMENSION_SLICE_PLUS_2_

#define _INDEX_HZ_OUTPUT_ _INDEX_OUTPUT_
#define _INDEX_HZ_OUTPUT_1D_(index) _INDEX_OUTPUT_1D_(index)

#else

! ---------------------------------------------------------------------------------
! Functions operating on horizontal slices of the domain only will not be vectorized
! ---------------------------------------------------------------------------------

#define _ARG_LOCATION_VARS_HZ_ _ARG_LOCATION_
#define _DECLARE_LOCATION_ARG_HZ_ _DECLARE_LOCATION_ARG_

! Spatial loop for quantities defined on horizontal slice of the full spatial domain.
#define _HORIZONTAL_LOOP_BEGIN_EX_(environment)
#define _HORIZONTAL_LOOP_END_

#define _DIMENSION_SLICE_HORIZONTAL_
#define _DIMENSION_SLICE_HORIZONTAL_PLUS_1_ ,dimension(:)
#define _DIMENSION_SLICE_HORIZONTAL_PLUS_2_ ,dimension(:,:)

#define _INDEX_HZ_OUTPUT_
#define _INDEX_HZ_OUTPUT_1D_(index) (index)

#endif

#define _HORIZONTAL_LOOP_BEGIN_ _HORIZONTAL_LOOP_BEGIN_EX_(environment)

! Expressions for indexing space-dependent FABM variables defined on the full spatial domain.
! These may be overridden by the host-specific driver (if it needs another order of dimensions).
! In that case, do not redefine the expressions here.
#ifndef _INDEX_ODE_
#define _INDEX_ODE_(variable) _INDEX_OUTPUT_1D_(variable)
#endif
#ifndef _INDEX_PPDD_
#define _INDEX_PPDD_(variable1,variable2) _INDEX_OUTPUT_2D_(variable1,variable2)
#endif
#ifndef _INDEX_CONSERVED_QUANTITY_
#define _INDEX_CONSERVED_QUANTITY_(variable) _INDEX_OUTPUT_1D_(variable)
#endif
#ifndef _INDEX_VERTICAL_MOVEMENT_
#define _INDEX_VERTICAL_MOVEMENT_(variable) _INDEX_OUTPUT_1D_(variable)
#endif
#ifndef _INDEX_SURFACE_FLUX_
#define _INDEX_SURFACE_FLUX_(index) _INDEX_HZ_OUTPUT_1D_(index)
#endif
#ifndef _INDEX_BOTTOM_FLUX_
#define _INDEX_BOTTOM_FLUX_(index) _INDEX_HZ_OUTPUT_1D_(index)
#endif

! For FABM: standard arguments used in calling biogeochemical routines.
#define _ARGUMENTS_ND_IN_ self%environment _ARG_LOCATION_ND_
#define _ARGUMENTS_IN_HZ_ self%environment _ARG_LOCATION_VARS_HZ_

! For BGC models: FABM arguments to routines implemented by biogeochemical models.
#define _ARGUMENTS_LOCAL_ environment _ARG_LOCATION_
#define _ARGUMENTS_ND_ environment _ARG_LOCATION_ND_
#define _ARGUMENTS_HZ_ environment _ARG_LOCATION_VARS_HZ_
#define _ARGUMENTS_DO_ _ARGUMENTS_ND_,rhs
#define _ARGUMENTS_DO_PPDD_ _ARGUMENTS_ND_,pp,dd
#define _ARGUMENTS_DO_SURFACE_ _ARGUMENTS_HZ_,flux
#define _ARGUMENTS_DO_BOTTOM_ _ARGUMENTS_HZ_,flux_pel,flux_ben
#define _ARGUMENTS_DO_BOTTOM_PPDD_ _ARGUMENTS_HZ_,pp,dd,benthos_offset
#define _ARGUMENTS_GET_VERTICAL_MOVEMENT_ _ARGUMENTS_ND_,velocity
#define _ARGUMENTS_GET_EXTINCTION_ _ARGUMENTS_ND_,extinction
#define _ARGUMENTS_GET_DRAG_ _ARGUMENTS_HZ_,drag
#define _ARGUMENTS_GET_ALBEDO_ _ARGUMENTS_HZ_,albedo
#define _ARGUMENTS_GET_CONSERVED_QUANTITIES_ _ARGUMENTS_ND_,sums
#define _ARGUMENTS_GET_HORIZONTAL_CONSERVED_QUANTITIES_ _ARGUMENTS_HZ_,sums
#define _ARGUMENTS_CHECK_STATE_ _ARGUMENTS_ND_,repair,valid
#define _ARGUMENTS_CHECK_SURFACE_STATE_ _ARGUMENTS_HZ_,repair,valid
#define _ARGUMENTS_CHECK_BOTTOM_STATE_ _ARGUMENTS_HZ_,repair,valid
#define _ARGUMENTS_INITIALIZE_STATE_ _ARGUMENTS_ND_
#define _ARGUMENTS_INITIALIZE_HORIZONTAL_STATE_ _ARGUMENTS_HZ_

! For BGC models: Declaration of FABM arguments to routines implemented by biogeochemical models.
#define _DECLARE_ARGUMENTS_LOCAL_ type (type_environment),intent(inout) :: environment;_DECLARE_LOCATION_ARG_
#define _DECLARE_ARGUMENTS_ND_ type (type_environment),intent(inout) :: environment;_DECLARE_LOCATION_ARG_ND_
#define _DECLARE_ARGUMENTS_HZ_ type (type_environment),intent(inout) :: environment;_DECLARE_LOCATION_ARG_HZ_
#define _DECLARE_ARGUMENTS_DO_  _DECLARE_ARGUMENTS_ND_;real(rk) _DIMENSION_SLICE_PLUS_1_,intent(inout) :: rhs
#define _DECLARE_ARGUMENTS_DO_PPDD_ _DECLARE_ARGUMENTS_ND_;real(rk) _DIMENSION_SLICE_PLUS_2_,intent(inout) :: pp,dd
#define _DECLARE_ARGUMENTS_DO_BOTTOM_ _DECLARE_ARGUMENTS_HZ_;real(rk) _DIMENSION_SLICE_HORIZONTAL_PLUS_1_,intent(inout) :: flux_pel,flux_ben
#define _DECLARE_ARGUMENTS_DO_BOTTOM_PPDD_ _DECLARE_ARGUMENTS_HZ_;real(rk) _DIMENSION_SLICE_HORIZONTAL_PLUS_2_,intent(inout) :: pp,dd;integer,intent(in) :: benthos_offset
#define _DECLARE_ARGUMENTS_DO_SURFACE_ _DECLARE_ARGUMENTS_HZ_;real(rk) _DIMENSION_SLICE_HORIZONTAL_PLUS_1_,intent(inout) :: flux
#define _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_ _DECLARE_ARGUMENTS_ND_;real(rk) _DIMENSION_SLICE_PLUS_1_,intent(inout) :: velocity
#define _DECLARE_ARGUMENTS_GET_EXTINCTION_ _DECLARE_ARGUMENTS_ND_;real(rk) _DIMENSION_SLICE_,intent(inout) :: extinction
#define _DECLARE_ARGUMENTS_GET_DRAG_ _DECLARE_ARGUMENTS_HZ_;real(rk) _DIMENSION_SLICE_HORIZONTAL_,intent(inout) :: drag
#define _DECLARE_ARGUMENTS_GET_ALBEDO_ _DECLARE_ARGUMENTS_HZ_;real(rk) _DIMENSION_SLICE_HORIZONTAL_,intent(inout) :: albedo
#define _DECLARE_ARGUMENTS_GET_CONSERVED_QUANTITIES_ _DECLARE_ARGUMENTS_ND_;real(rk) _DIMENSION_SLICE_PLUS_1_,intent(inout) :: sums
#define _DECLARE_ARGUMENTS_GET_HORIZONTAL_CONSERVED_QUANTITIES_ _DECLARE_ARGUMENTS_HZ_;real(rk) _DIMENSION_SLICE_HORIZONTAL_PLUS_1_,intent(inout) :: sums
#define _DECLARE_ARGUMENTS_CHECK_STATE_ _DECLARE_ARGUMENTS_ND_;logical,intent(in) :: repair;logical,intent(inout) :: valid
#define _DECLARE_ARGUMENTS_CHECK_SURFACE_STATE_ _DECLARE_ARGUMENTS_HZ_;logical,intent(in) :: repair;logical,intent(inout) :: valid
#define _DECLARE_ARGUMENTS_CHECK_BOTTOM_STATE_ _DECLARE_ARGUMENTS_HZ_;logical,intent(in) :: repair;logical,intent(inout) :: valid
#define _DECLARE_ARGUMENTS_INITIALIZE_STATE_ _DECLARE_ARGUMENTS_ND_
#define _DECLARE_ARGUMENTS_INITIALIZE_HORIZONTAL_STATE_ _DECLARE_ARGUMENTS_HZ_

! For backward compatibility (pre 20 June 2013)
#define _FABM_ARGS_DO_RHS_ _ARGUMENTS_DO_
#define _FABM_ARGS_DO_PPDD_ _ARGUMENTS_DO_PPDD_
#define _FABM_ARGS_DO_BENTHOS_RHS_ _ARGUMENTS_DO_BOTTOM_
#define _FABM_ARGS_DO_BENTHOS_PPDD_ _ARGUMENTS_DO_BOTTOM_PPDD_
#define _FABM_ARGS_GET_SURFACE_EXCHANGE_ _ARGUMENTS_DO_SURFACE_
#define _FABM_ARGS_GET_EXTINCTION_ _ARGUMENTS_GET_EXTINCTION_
#define _FABM_ARGS_GET_DRAG_ _ARGUMENTS_GET_DRAG_
#define _FABM_ARGS_GET_ALBEDO_ _ARGUMENTS_GET_ALBEDO_
#define _FABM_ARGS_GET_VERTICAL_MOVEMENT_ _ARGUMENTS_GET_VERTICAL_MOVEMENT_
#define _FABM_ARGS_GET_CONSERVED_QUANTITIES_ _ARGUMENTS_GET_CONSERVED_QUANTITIES_
#define _FABM_ARGS_CHECK_STATE_ _ARGUMENTS_CHECK_STATE_

! For backward compatibility (pre 20 June 2013)
#define _DECLARE_FABM_ARGS_DO_RHS_  _DECLARE_ARGUMENTS_DO_
#define _DECLARE_FABM_ARGS_DO_PPDD_ _DECLARE_ARGUMENTS_DO_PPDD_
#define _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_ _DECLARE_ARGUMENTS_DO_SURFACE_
#define _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_ _DECLARE_ARGUMENTS_DO_BOTTOM_
#define _DECLARE_FABM_ARGS_DO_BENTHOS_PPDD_ _DECLARE_ARGUMENTS_DO_BOTTOM_PPDD_
#define _DECLARE_FABM_ARGS_GET_VERTICAL_MOVEMENT_ _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
#define _DECLARE_FABM_ARGS_GET_EXTINCTION_ _DECLARE_ARGUMENTS_GET_EXTINCTION_
#define _DECLARE_FABM_ARGS_GET_DRAG_ _DECLARE_ARGUMENTS_GET_DRAG_
#define _DECLARE_FABM_ARGS_GET_ALBEDO_ _DECLARE_ARGUMENTS_GET_ALBEDO_
#define _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_ _DECLARE_ARGUMENTS_GET_CONSERVED_QUANTITIES_
#define _DECLARE_FABM_ARGS_CHECK_STATE_ _DECLARE_ARGUMENTS_CHECK_STATE_

! Macros for declaring/accessing variable identifiers of arbitrary type.
#define _TYPE_STATE_VARIABLE_ID_ type (type_state_variable_id)
#define _TYPE_DIAGNOSTIC_VARIABLE_ID_ type (type_diagnostic_variable_id)
#define _TYPE_DEPENDENCY_ID_ type (type_dependency_id)
#define _TYPE_CONSERVED_QUANTITY_ID_ type (type_conserved_quantity_id)

! For BGC models: Expressions for setting space-dependent FABM variables defined on the full spatial domain.
#define _SET_ODE_(variable,value) rhs _INDEX_ODE_(variable%state_index) = rhs _INDEX_ODE_(variable%state_index) + (value)/self%dt
#define _SET_ODE_BEN_(variable,value) flux_ben _INDEX_BOTTOM_FLUX_(variable%bottom_state_index) = flux_ben _INDEX_BOTTOM_FLUX_(variable%bottom_state_index) + (value)/self%dt
#define _SET_BOTTOM_EXCHANGE_(variable,value) flux_pel _INDEX_BOTTOM_FLUX_(variable%state_index) = flux_pel _INDEX_BOTTOM_FLUX_(variable%state_index) + (value)/self%dt
#define _SET_SURFACE_EXCHANGE_(variable,value) flux _INDEX_SURFACE_FLUX_(variable%state_index) = value/self%dt
#define _SET_ODE_(variable,value) rhs _INDEX_ODE_(variable%state_index) = rhs _INDEX_ODE_(variable%state_index) + (value)/self%dt
#define _SET_DD_(variable1,variable2,value) dd _INDEX_PPDD_(variable1%state_index,variable2%state_index) = dd _INDEX_PPDD_(variable1%state_index,variable2%state_index) + (value)/self%dt
#define _SET_PP_(variable1,variable2,value) pp _INDEX_PPDD_(variable1%state_index,variable2%state_index) = pp _INDEX_PPDD_(variable1%state_index,variable2%state_index) + (value)/self%dt
#define _SET_EXTINCTION_(value) extinction _INDEX_OUTPUT_ = extinction _INDEX_OUTPUT_ + (value)
#define _SCALE_DRAG_(value) drag _INDEX_HZ_OUTPUT_ = drag _INDEX_HZ_OUTPUT_ * (value)
#define _SET_ALBEDO_(value) albedo _INDEX_HZ_OUTPUT_ = albedo _INDEX_HZ_OUTPUT_ + (value)
#define _SET_CONSERVED_QUANTITY_(variable,value) sums _INDEX_CONSERVED_QUANTITY_(variable%cons_index) = sums _INDEX_CONSERVED_QUANTITY_(variable%cons_index) + (value)
#define _SET_VERTICAL_MOVEMENT_(variable,value) velocity _INDEX_VERTICAL_MOVEMENT_(variable%state_index) = value/self%dt
#define _INVALIDATE_STATE_ valid = .false.
#define _REPAIR_STATE_ repair

! For BGC models: quick expressions for setting a single element in both the destruction and production matrix.
#define _SET_DD_SYM_(variable1,variable2,value) _SET_DD_(variable1,variable2,value);_SET_PP_(variable2,variable1,value)
#define _SET_PP_SYM_(variable1,variable2,value) _SET_PP_(variable1,variable2,value);_SET_DD_(variable2,variable1,value)

! For BGC models: macro to determine whether a variable identifier is in use (i.e., has been registered with FABM)
#define _VARIABLE_REGISTERED_(variable) associated(variable%link)
#define _AVAILABLE_(variable) associated(variable%data%p)

! Within FABM: read/write variable access.
#define _GET_EX_(variable,target) target = variable%p _INDEX_LOCATION_
#define _GET_HORIZONTAL_EX_(variable,target) target = variable%p _INDEX_LOCATION_HZ_
#define _GET_GLOBAL_EX_(variable,target) target = variable%p
#define _SET_EX_(variable,value) variable%p _INDEX_LOCATION_ = value
#define _SET_HORIZONTAL_EX_(variable,value) variable%p _INDEX_LOCATION_HZ_ = value
#define _SET_GLOBAL_EX_(variable,value) variable%p = value

! For BGC models: read/write variable access.
#define _GET_(variable,target) _GET_EX_(variable%data,target)
#define _GET_HORIZONTAL_(variable,target) _GET_HORIZONTAL_EX_(variable%horizontal_data,target)
#define _GET_GLOBAL_(variable,target) _GET_GLOBAL_EX_(variable%global_data,target)
#define _SET_(variable,value) _SET_EX_(variable%data,value)
#define _SET_HORIZONTAL_(variable,value) _SET_HORIZONTAL_EX_(variable%horizontal_data,value)
#define _SET_GLOBAL_(variable,value) _SET_GLOBAL_EX_(variable%global_data,value)
#define _SET_DIAGNOSTIC_(variable,value) environment%diag(_PREARG_LOCATION_ variable%diag_index) = value
#define _SET_HORIZONTAL_DIAGNOSTIC_(variable,value) environment%diag_hz(_PREARG_LOCATION_HZ_ variable%horizontal_diag_index) = value

! For backward compatibility: old macros to access variable data.
#define _GET_DEPENDENCY_(variable,target) _GET_(variable,target)
#define _GET_DEPENDENCY_HZ_(variable,target) _GET_HORIZONTAL_(variable,target)
#define _GET_DEPENDENCY_SCALAR_(variable,target) _GET_GLOBAL_(variable,target)
#define _GET_STATE_(variable,target) _GET_(variable,target)
#define _GET_STATE_BEN_(variable,target) _GET_HORIZONTAL_(variable,target)
#define _SET_STATE_(variable,target) _SET_(variable,target)
#define _SET_STATE_BEN_(variable,target) _SET_HORIZONTAL_(variable,target)
#define _GET_STATE_EX_(env,variable,target) _GET_EX_(variable,target)
#define _GET_STATE_BEN_EX_(env,variable,target) _GET_HORIZONTAL_EX_(variable,target)
#define _SET_STATE_EX_(env,variable,value) _SET_EX_(variable,value)
#define _SET_STATE_BEN_EX_(env,variable,value) _SET_HORIZONTAL_EX_(variable,value)
#define _SET_DIAG_(variable,value) _SET_DIAGNOSTIC_(variable,value)
#define _SET_DIAG_HZ_(variable,value) _SET_HORIZONTAL_DIAGNOSTIC_(variable,value)
#define _FABM_HZ_LOOP_BEGIN_ _HORIZONTAL_LOOP_BEGIN_
#define _FABM_HZ_LOOP_END_ _HORIZONTAL_LOOP_END_

! Work-in-progress: extra definitions for coupling to pure-1D models [ERSEM]
! Currently these are GOTM-specific - more logic will be needed to set these to
! appropriate values for non-column or non-vectorized models.
#define _DOMAIN_1D_ fabm_loop_start:fabm_loop_stop
#define _GET_STATE_1D_(variable,target) target = environment%var(variable%dependencyid)%data(_DOMAIN_1D_)
#define _GET_DEPENDENCY_1D_(variable,target) target = environment%var(variable)%data(_DOMAIN_1D_)
#define _LOOP_BEGIN_1D_
#define _LOOP_END_1D_
#define _HORIZONTAL_LOOP_BEGIN_1D_
#define _HORIZONTAL_LOOP_END_1D_
#ifndef _INDEX_ODE_1D_
#define _INDEX_ODE_1D_(variable) (1:fabm_loop_stop-fabm_loop_start+1,variable)
#endif
#define _SET_ODE_1D_(variable,value) rhs _INDEX_ODE_1D_(variable%id) = rhs _INDEX_ODE_1D_(variable%id) + (value)
#define _SET_EXTINCTION_1D_(value) extinction(1:fabm_loop_stop-fabm_loop_start+1) = extinction(fabm_loop_start:fabm_loop_stop) + value
#define _SET_VERTICAL_MOVEMENT_1D_(variable,value) velocity(1:fabm_loop_stop-fabm_loop_start+1,variable%id) = value

! For the definitions below, it is assumed that the vertical dimension is vectorized!
#define _DOMAIN_HZ_1D_ 1
#define _INDEX_HZ_1D_ _VARIABLE_1DLOOP_
#define _GET_STATE_HZ_1D_(variable,target) target = environment%var(variable%dependencyid)%data(_INDEX_HZ_1D_)
#define _GET_DEPENDENCY_HZ_1D_(variable,target) target = environment%var(variable)%data(_INDEX_HZ_1D_)
#define _GET_STATE_BEN_1D_(variable,target) target = environment%var_hz(variable%dependencyid)%data
#define _SET_BOTTOM_FLUX_1D_(variable,value) flux_pel(variable%id) = flux_pel(variable%id) + (value)
#define _SET_ODE_BEN_1D_(variable,value) flux_ben(variable%id) = flux_ben(variable%id) + (value)

! For backward compatibility (pre 20 June 2013)
#define _FABM_LOOP_BEGIN_ _LOOP_BEGIN_
#define _FABM_LOOP_END_ _LOOP_END_
#define _FABM_HORIZONTAL_LOOP_BEGIN_ _HORIZONTAL_LOOP_BEGIN_
#define _FABM_HORIZONTAL_LOOP_END_ _HORIZONTAL_LOOP_END_

#define _GET_WITHOUT_BACKGROUND_(variable,target) target = max(0.0_rk,variable%data%p _INDEX_LOCATION_-variable%background)

#ifdef _FABM_DEPTH_DIMENSION_INDEX_

! ---------------------------------------------------------------------------------
! The model has a vertical dimension.
! Functions operating over all values in the vertical will be vectorized.
! ---------------------------------------------------------------------------------

! Identify variable that contains index of vertical dimension
#if _FABM_DEPTH_DIMENSION_INDEX_==1
#define _VARIABLE_VERT_LOOP_ i__
#elif _FABM_DEPTH_DIMENSION_INDEX_==2
#define _VARIABLE_VERT_LOOP_ j__
#else
#define _VARIABLE_VERT_LOOP_ k__
#endif

! Define arguments that describe the horizontal location (fixed during vertical loops)
#if _FABM_DIMENSION_COUNT_==1
#define _ARG_LOCATION_VERTICAL_LOOP_
#else
#define _ARG_LOCATION_VERTICAL_LOOP_ ,_LOCATION_HZ_
#endif

! Define arguments for routines operating on the vertical dimension.
#define _ARG_LOCATION_VERT_ ,fabm_loop_start,fabm_loop_stop _ARG_LOCATION_VERTICAL_LOOP_
#define _DECLARE_LOCATION_ARG_VERT_ integer,intent(in) :: fabm_loop_start,fabm_loop_stop _ARG_LOCATION_VERTICAL_LOOP_;integer :: _VARIABLE_VERT_LOOP_

! Define statements to start and finish vertical loops
#ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
#define _VERTICAL_LOOP_BEGIN_ do _VARIABLE_VERT_LOOP_=fabm_loop_stop,fabm_loop_start,-1
#else
#define _VERTICAL_LOOP_BEGIN_ do _VARIABLE_VERT_LOOP_=fabm_loop_start,fabm_loop_stop
#endif
#define _VERTICAL_LOOP_END_ end do

#else

! ---------------------------------------------------------------------------------
! The model lacks a vertical dimension.
! Functions operating over all values in the vertical will not be vectorized.
! ---------------------------------------------------------------------------------

#define _ARG_LOCATION_VERT_ _ARG_LOCATION_
#define _DECLARE_LOCATION_ARG_VERT_ _DECLARE_LOCATION_ARG_

#define _VERTICAL_LOOP_BEGIN_
#define _VERTICAL_LOOP_END_

#endif

#define _ARGUMENTS_VERT_ environment _ARG_LOCATION_VERT_
#define _DECLARE_ARGUMENTS_VERT_ type (type_environment),intent(inout) :: environment;_DECLARE_LOCATION_ARG_VERT_

! Make sure that the primary preprocessor macro has been defined.
#ifndef _FABM_DIMENSION_COUNT_
#error Preprocessor variable _FABM_DIMENSION_COUNT_ must be defined.
#endif

! If the spatial context is 0D, some preprocessor variables can be inferred.
#if _FABM_DIMENSION_COUNT_==0
#define _LOCATION_
#define _LOCATION_DIMENSIONS_
#define _FABM_HORIZONTAL_IS_SCALAR_
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

! Constants related to floating point precision; used throughout FABM.
#define _NP_ 8
#define REALTYPE real(kind=_NP_)
#define _ZERO_ 0.0__NP_
#define _ONE_  1.0__NP_

! Older Fortran compilers do not allow derived types to contain allocatable members
! (A Fortran >95 feature, defined in ISO Technical Report TR 15581 and part of the Fortran 2003 specification).
! As a workaround, they can be declared with the pointer attribute, which does bring a slight performance penalty.
! By using the below preprocessor definitions, the allocatable attribute is automatically replaced by the pointer
! attribute where needed, and related function calls are changed as well.
! These are used through FABM.
#ifdef _ISO_TR_15581_
#define _ALLOCATABLE_ allocatable
#define _NULL_ 
#define _ALLOCATED_ allocated
#else
#define _ALLOCATABLE_ pointer
#define _NULL_ =>null()
#define _ALLOCATED_ associated
#endif

! Data type for location variable(s); only used within this file.
#define _LOCATION_TYPE_ integer

! Define dimension attribute and index specifyer for horizontal (2D) fields.
#ifdef _FABM_HORIZONTAL_IS_SCALAR_
#define _INDEX_LOCATION_HZ_
#define _ATTR_LOCATION_DIMENSIONS_HZ_
#define _ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE_ ,dimension(:)
#define _ARG_LOCATION_HZ_
#define _ARG_LOCATION_DIMENSIONS_HZ_
#else
#define _INDEX_LOCATION_HZ_ (_LOCATION_HZ_)
#define _ATTR_LOCATION_DIMENSIONS_HZ_ ,dimension(_LOCATION_DIMENSIONS_HZ_)
#define _ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE_ ,dimension(:,_LOCATION_DIMENSIONS_HZ_)
#define _ARG_LOCATION_HZ_ ,_LOCATION_HZ_
#define _ARG_LOCATION_DIMENSIONS_HZ_ ,_LOCATION_DIMENSIONS_HZ_
#endif

! Define dimension attribute and index specifyer for full 3D fields.
#if _FABM_DIMENSION_COUNT_>0
#define _INDEX_LOCATION_ (_LOCATION_)
#define _ATTR_LOCATION_DIMENSIONS_ ,dimension(_LOCATION_DIMENSIONS_)
#define _ATTR_LOCATION_DIMENSIONS_PLUS_ONE_ ,dimension(:,_LOCATION_DIMENSIONS_)
#define _ARG_LOCATION_ ,_LOCATION_
#define _ARG_LOCATION_DIMENSIONS_ ,:
#define _DECLARE_LOCATION_ARG_ _LOCATION_TYPE_,intent(in) :: _LOCATION_
#else
#define _INDEX_LOCATION_
#define _ATTR_LOCATION_DIMENSIONS_
#define _ATTR_LOCATION_DIMENSIONS_PLUS_ONE_ ,dimension(:)
#define _ARG_LOCATION_
#define _ARG_LOCATION_DIMENSIONS_
#define _DECLARE_LOCATION_ARG_
#endif

#ifdef _FABM_USE_1D_LOOP_

! 1D vectorized: FABM subroutines operate on one spatial dimension.

! Make sure there is at least one dimension to vectorize.
#if _FABM_DIMENSION_COUNT_<1
#error Cannot build 1D vectorized version of FABM (_FABM_USE_1D_LOOP_) with the total number of spatial dimensions (_FABM_DIMENSION_COUNT_) being less than 1.
#endif

! For a 1D spatial domain, the vectorzied version of FABM must iterate over that sole dimension.
! Therefore, the _LOCATION_1DLOOP_ and _VARIABLE_1DLOOP_ macros can be inferred.
#if _FABM_DIMENSION_COUNT_>1
#define _ARG_LOCATION_1DLOOP_ ,_LOCATION_1DLOOP_
#else
#define _ARG_LOCATION_1DLOOP_
#define _LOCATION_1DLOOP_
#define _VARIABLE_1DLOOP_ _LOCATION_
#endif

! Check for additional preprocessor macros that a required when _FABM_USE_1D_LOOP_ is defined.
#ifndef _VARIABLE_1DLOOP_
#error Building 1D vectorized version of FABM: preprocessor variable _VARIABLE_1DLOOP_ must be defined.
#endif
#ifndef _LOCATION_1DLOOP_
#error Building 1D vectorized version of FABM: preprocessor variable _LOCATION_1DLOOP_ must be defined.
#endif

! Dummy argument and argument declaration for location specification.
#define _ARG_LOCATION_ND_ ,fabm_loop_start,fabm_loop_stop _ARG_LOCATION_1DLOOP_
#define _DECLARE_LOCATION_ARG_ND_ _LOCATION_TYPE_,intent(in) :: fabm_loop_start,fabm_loop_stop _ARG_LOCATION_1DLOOP_; _LOCATION_TYPE_ :: _VARIABLE_1DLOOP_

! Beginning and end of spatial loop
#define _FABM_LOOP_BEGIN_ do _VARIABLE_1DLOOP_=fabm_loop_start,fabm_loop_stop
#define _FABM_LOOP_END_ end do

! Dimensionality of generic space-dependent arguments.
#define _ATTR_DIMENSIONS_0_ ,dimension(:)
#define _ATTR_DIMENSIONS_1_ ,dimension(:,:)
#define _ATTR_DIMENSIONS_2_ ,dimension(:,:,:)

! Expressions for indexing space-dependent FABM variables defined on the full spatial domain.
! These may be overridden by the host-specific driver (if it needs another order of dimensions).
! In that case, do not redefine the expressions here.
#ifndef _INDEX_ODE_
#define _INDEX_ODE_(variable) (_VARIABLE_1DLOOP_-fabm_loop_start+1,variable)
#endif
#ifndef _INDEX_PPDD_
#define _INDEX_PPDD_(variable1,variable2) (_VARIABLE_1DLOOP_-fabm_loop_start+1,variable1,variable2)
#endif
#ifndef _INDEX_CONSERVED_QUANTITY_
#define _INDEX_CONSERVED_QUANTITY_(variable) (_VARIABLE_1DLOOP_-fabm_loop_start+1,variable)
#endif
#ifndef _INDEX_VERTICAL_MOVEMENT_
#define _INDEX_VERTICAL_MOVEMENT_(variable) (_VARIABLE_1DLOOP_-fabm_loop_start+1,variable)
#endif
#define _INDEX_EXTINCTION_ (_VARIABLE_1DLOOP_-fabm_loop_start+1)

#else

! Not vectorized: FABM subroutines operate one the local state only.

! Dummy argument and argument declaration for location specification.
#define _ARG_LOCATION_ND_ _ARG_LOCATION_
#define _DECLARE_LOCATION_ARG_ND_ _DECLARE_LOCATION_ARG_

! Beginning and end of spatial loop
#define _FABM_LOOP_BEGIN_
#define _FABM_LOOP_END_

! Dimensionality of generic space-dependent arguments.
#define _ATTR_DIMENSIONS_0_
#define _ATTR_DIMENSIONS_1_ ,dimension(:)
#define _ATTR_DIMENSIONS_2_ ,dimension(:,:)

! Expressions for indexing space-dependent FABM variables defined on the full spatial domain.
#define _INDEX_ODE_(variable) (variable)
#define _INDEX_PPDD_(variable1,variable2) (variable1,variable2)
#define _INDEX_EXTINCTION_
#define _INDEX_CONSERVED_QUANTITY_(variable) (variable)
#define _INDEX_VERTICAL_MOVEMENT_(variable) (variable)

#endif

#if defined(_FABM_USE_1D_LOOP_)&&(!defined(_FABM_1D_LOOP_IN_VERTICAL_))

! Functions operating on horizontonal slices of the domain will be vectorized just like the full domain.

#define _ARG_LOCATION_VARS_HZ_ _ARG_LOCATION_ND_
#define _DECLARE_LOCATION_ARG_HZ_ _DECLARE_LOCATION_ARG_ND_

! Spatial loop for quantities defined on horizontal slice of the full spatial domain.
#define _FABM_HZ_LOOP_BEGIN_ _FABM_LOOP_BEGIN_
#define _FABM_HZ_LOOP_END_ _FABM_LOOP_END_

! Vertical dimension is not among those vectorized:
! dimensionality of horizontal arrays will be equal to that of full domain arrays.
#define _ATTR_DIMENSIONS_0_HZ_ _ATTR_DIMENSIONS_0_
#define _ATTR_DIMENSIONS_1_HZ_ _ATTR_DIMENSIONS_1_
#define _ATTR_DIMENSIONS_2_HZ_ _ATTR_DIMENSIONS_2_

! Expressions for indexing space-dependent FABM variables defined on horizontal slices of the domain.
! May be overridden by host to reverse the dimension order.
#ifndef _INDEX_SURFACE_EXCHANGE_
#define _INDEX_SURFACE_EXCHANGE_(index) (_VARIABLE_1DLOOP_-fabm_loop_start+1,index)
#endif

#else

! Functions operating on horizontal slices of the domain only will not be vectorized

#define _ARG_LOCATION_VARS_HZ_ _ARG_LOCATION_
#define _DECLARE_LOCATION_ARG_HZ_ _DECLARE_LOCATION_ARG_

! Spatial loop for quantities defined on horizontal slice of the full spatial domain.
#define _FABM_HZ_LOOP_BEGIN_
#define _FABM_HZ_LOOP_END_

! Expressions for indexing space-dependent FABM variables defined on horizontal slices of the domain.
#define _INDEX_SURFACE_EXCHANGE_(index) (index)

#define _ATTR_DIMENSIONS_0_HZ_
#define _ATTR_DIMENSIONS_1_HZ_ ,dimension(:)
#define _ATTR_DIMENSIONS_2_HZ_ ,dimension(:,:)

#endif


! For FABM: standard arguments used in calling biogeochemical routines.
#define _FABM_ARGS_ND_IN_ root%environment _ARG_LOCATION_ND_
#define _FABM_ARGS_IN_HZ_ root%environment _ARG_LOCATION_VARS_HZ_

! For BGC models: FABM arguments to routines implemented by biogeochemical models.
#define _FABM_ARGS_ND_ environment _ARG_LOCATION_ND_
#define _FABM_ARGS_HZ_ environment _ARG_LOCATION_VARS_HZ_
#define _FABM_ARGS_DO_RHS_ _FABM_ARGS_ND_,rhs
#define _FABM_ARGS_DO_PPDD_ _FABM_ARGS_ND_,pp,dd
#define _FABM_ARGS_GET_EXTINCTION_ _FABM_ARGS_ND_,extinction
#define _FABM_ARGS_GET_CONSERVED_QUANTITIES_ _FABM_ARGS_ND_,sums
#define _FABM_ARGS_GET_SURFACE_EXCHANGE_ _FABM_ARGS_HZ_,flux
#define _FABM_ARGS_CHECK_STATE_ _FABM_ARGS_ND_,repair,valid

! For BGC models: Declaration of FABM arguments to routines implemented by biogeochemical models.
#define _DECLARE_FABM_ARGS_ND_ type (type_environment),intent(inout) :: environment;_DECLARE_LOCATION_ARG_ND_
#define _DECLARE_FABM_ARGS_HZ_ type (type_environment),intent(inout) :: environment;_DECLARE_LOCATION_ARG_HZ_
#define _DECLARE_FABM_ARGS_DO_RHS_  _DECLARE_FABM_ARGS_ND_;REALTYPE _ATTR_DIMENSIONS_1_,intent(inout) :: rhs
#define _DECLARE_FABM_ARGS_DO_PPDD_ _DECLARE_FABM_ARGS_ND_;REALTYPE _ATTR_DIMENSIONS_2_,intent(inout) :: pp,dd
#define _DECLARE_FABM_ARGS_GET_EXTINCTION_ _DECLARE_FABM_ARGS_ND_;REALTYPE _ATTR_DIMENSIONS_0_,intent(inout) :: extinction
#define _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_ _DECLARE_FABM_ARGS_ND_;REALTYPE _ATTR_DIMENSIONS_1_,intent(inout) :: sums
#define _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_ _DECLARE_FABM_ARGS_HZ_;REALTYPE _ATTR_DIMENSIONS_1_HZ_,intent(inout) :: flux
#define _DECLARE_FABM_ARGS_CHECK_STATE_ _DECLARE_FABM_ARGS_ND_;logical,intent(in) :: repair;logical,intent(inout) :: valid

! Macros for declaring/accessing variable identifiers of arbitrary type.
#define _TYPE_STATE_VARIABLE_ID_ type (type_state_variable_id)
#define _TYPE_DIAGNOSTIC_VARIABLE_ID_ integer
#define _TYPE_DEPENDENCY_ID_ integer
#define _TYPE_CONSERVED_QUANTITY_ID_ integer

! For BGC models: Expressions for setting space-dependent FABM variables defined on the full spatial domain.
#define _SET_ODE_(variable,value) rhs _INDEX_ODE_(variable%id) = rhs _INDEX_ODE_(variable%id) + (value)
#define _SET_DD_(variable1,variable2,value) dd _INDEX_PPDD_(variable1%id,variable2%id) = dd _INDEX_PPDD_(variable1%id,variable2%id) + (value)
#define _SET_PP_(variable1,variable2,value) pp _INDEX_PPDD_(variable1%id,variable2%id) = pp _INDEX_PPDD_(variable1%id,variable2%id) + (value)
#define _SET_EXTINCTION_(value) extinction _INDEX_EXTINCTION_ = extinction _INDEX_EXTINCTION_ + (value)
#define _SET_CONSERVED_QUANTITY_(variable,value) sums _INDEX_CONSERVED_QUANTITY_(variable) = sums _INDEX_CONSERVED_QUANTITY_(variable) + (value)
#define _SET_VERTICAL_MOVEMENT_(variable,value) vertical_movement _INDEX_VERTICAL_MOVEMENT_(variable%id) = value
#define _SET_SURFACE_EXCHANGE_(variable,value) flux _INDEX_SURFACE_EXCHANGE_(variable%id) = value
#define _INVALIDATE_STATE_ valid = .false.
#define _REPAIR_STATE_ repair

! For BGC models: quick expressions for setting a single element in both the destruction and production matrix.
#define _SET_DD_SYM_(variable1,variable2,value) _SET_DD_(variable1,variable2,value);_SET_PP_(variable1,variable2,value)
#define _SET_PP_SYM_(variable1,variable2,value) _SET_PP_(variable1,variable2,value);_SET_DD_(variable1,variable2,value)

! For BGC models: read-only access to values of external dependencies
#define _GET_DEPENDENCY_(variable,target) target = environment%var(variable)%data _INDEX_LOCATION_
#define _GET_DEPENDENCY_HZ_(variable,target) target = environment%var_hz(variable)%data _INDEX_LOCATION_HZ_

! For FABM: read/write access to state variables
#define _GET_STATE_EX_(env,variable,target) target = env%var(variable%dependencyid)%data _INDEX_LOCATION_
#define _SET_STATE_EX_(env,variable,value) env%var(variable%dependencyid)%data _INDEX_LOCATION_ = value

! For BGC models: read/write access to state variable values
#define _GET_STATE_(variable,target) _GET_STATE_EX_(environment,variable,target)
#define _SET_STATE_(variable,target) _SET_STATE_EX_(environment,variable,target)

! For BGC models: write access to diagnostic variables
#ifdef _FABM_MANAGE_DIAGNOSTICS_
#define _SET_DIAG_(index,value) if (index.ne.id_not_used) environment%diag(index,_LOCATION_) = value
#define _SET_DIAG_HZ_(index,value) if (index.ne.id_not_used) environment%diag_hz(index) = value
#else
#define _SET_DIAG_(index,value) if (index.ne.id_not_used) environment%var(index)%data _INDEX_LOCATION_ = value
#define _SET_DIAG_HZ_(index,value) if (index.ne.id_not_used) environment%var_hz(index)%data _INDEX_LOCATION_HZ_ = value
#endif


! Work-in-progress: extra definitions for coupling to pure-1D models [ERSEM]
#define _GET_NLEV_ 150
#define _GET_STATE_1D_(variable,target) target(fabm_loop_start:fabm_loop_stop) = env%var(variable%dependencyid)%data(fabm_loop_start:fabm_loop_stop)
#define _GET_DEPENDENCY_1D_(variable,target) target(fabm_loop_start:fabm_loop_stop) = env%var(variable)%data(fabm_loop_start:fabm_loop_stop)
#define _FABM_LOOP_BEGIN_1D_
#define _FABM_LOOP_END_1D_
#ifndef _INDEX_ODE_1D_
#define _INDEX_ODE_1D_(variable) (fabm_loop_start:fabm_loop_stop,variable)
#endif
#define _SET_ODE_1D_ rhs _INDEX_ODE_1D_(variable%id) = rhs _INDEX_ODE_1D_(variable%id) + (value)
#define _SET_EXTINCTION_1D_ extinction(fabm_loop_start:fabm_loop_stop) = extinction(fabm_loop_start:fabm_loop_stop) + value(fabm_loop_start:fabm_loop_stop)

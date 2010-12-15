#define REALTYPE double precision
#define _ZERO_ 0.0d0
#define _ONE_  1.0d0

! Preprocessor definition below have been taken from mom4p0d

#if defined(__INTEL_COMPILER) || defined(__IBMC__)
#define _F95
#define _F2000
#endif

#ifdef _F95
!DEC$ MESSAGE:'Using PURE'
#define _PURE PURE
#else
!DEC$ MESSAGE:'Not using PURE'
#define _PURE 
#endif

! Older Fortran compilers do not allow derived types to contain allocatable members (A Fortran >95 feature).
! As a workaround, they can be declared with the pointer attribute. By using the below preprocessor definitions,
! the allocatable attribute is automatically repalced by the pointer attribute where needed, and related function
! calls are chanegd as well.
#ifdef _F2000
!DEC$ MESSAGE:'Converting pointers to allocatable components'
#define _ALLOCATABLE ALLOCATABLE
#define _NULL 
#define _ALLOCATED ALLOCATED
#else
!DEC$ MESSAGE:'Using pointers'
#define _ALLOCATABLE POINTER
#define _NULL =>NULL()
#define _ALLOCATED ASSOCIATED
#endif

! Data type for location variable(s)
#define LOCATION_TYPE integer

! Define remaining dimensions in a 1D loop.
! NB if there is only one spatial dimension, there are no remaining dimensions!
#if RMBM_DIMENSIONS>1
#define ARG_LOCATION_1DLOOP ,LOCATION_1DLOOP
#else
#define ARG_LOCATION_1DLOOP
#define LOCATION_1DLOOP
#define VARIABLE_1DLOOP LOCATION
#endif

! Define dimension attribute and index specifyer for horizontal (2D) fields.
#ifdef RMBM_HORIZONTAL_IS_SCALAR
#define INDEX_LOCATION_HZ
#define ATTR_LOCATION_DIMENSIONS_HZ
#define ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE ,dimension(:)
#define ARG_LOCATION_HZ
#define ARG_LOCATION_DIMENSIONS_HZ
#else
#define INDEX_LOCATION_HZ (LOCATION_HZ)
#define ATTR_LOCATION_DIMENSIONS_HZ ,dimension(LOCATION_DIMENSIONS_HZ)
#define ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE ,dimension(:,LOCATION_DIMENSIONS_HZ)
#define ARG_LOCATION_HZ ,LOCATION_HZ
#define ARG_LOCATION_DIMENSIONS_HZ ,LOCATION_DIMENSIONS_HZ
#endif

! Define dimension attribute and index specifyer for full 3D fields.
#define INDEX_LOCATION (LOCATION)
#define ATTR_LOCATION_DIMENSIONS ,dimension(LOCATION_DIMENSIONS)
#define ATTR_LOCATION_DIMENSIONS_PLUS_ONE ,dimension(:,LOCATION_DIMENSIONS)
#define ARG_LOCATION ,LOCATION
#define ARG_LOCATION_DIMENSIONS ,:

#define RMBM_ARGS_0D environment ARG_LOCATION
#define DECLARE_RMBM_ARGS_0D type (type_environment),intent(inout) :: environment;LOCATION_TYPE,intent(in) :: LOCATION
#define RMBM_ARGS_IN_0D root%environment ARG_LOCATION

! Spatial loop for quantities defined on hortizontal slice of the full spatial domain.
#define _RMBM_ENTER_HZ_
#define _RMBM_LEAVE_HZ_

! Expressions for indexing space-dependent RMBM variables defined on horizontal slices of the domain.
#define _INDEX_SURFACE_EXCHANGE_(index) (index)

#ifdef RMBM_USE_1D_LOOP

! 1D vectorized: RMBM subroutines operate on one spatial dimension.

! Dummy argument and argument declaration for location specification.
#define LOCATION_ND rmbm_loop_start,rmbm_loop_stop ARG_LOCATION_1DLOOP
#define DECLARE_LOCATION_ARG_ND LOCATION_TYPE,intent(in) :: rmbm_loop_start,rmbm_loop_stop ARG_LOCATION_1DLOOP; LOCATION_TYPE :: VARIABLE_1DLOOP

! Beginning and end of spatial loop
#define _RMBM_LOOP_BEGIN_ do VARIABLE_1DLOOP=rmbm_loop_start,rmbm_loop_stop
#define _RMBM_LOOP_END_ end do

! Dimensionality of generic space-dependent arguments.
#define ATTR_DIMENSIONS_0 ,dimension(:)
#define ATTR_DIMENSIONS_1 ,dimension(:,:)
#define ATTR_DIMENSIONS_2 ,dimension(:,:,:)

! Expressions for indexing space-dependent RMBM variables defined on the full spatial domain.
! These may be overridden by the host-specific driver (if it needs another order of dimensions).
! In that case, do not redefine the expressions here.
#ifndef _INDEX_ODE_
#define _INDEX_ODE_(variable) (VARIABLE_1DLOOP,variable)
#endif
#ifndef _INDEX_PPDD_
#define _INDEX_PPDD_(variable1,variable2) (VARIABLE_1DLOOP,variable1,variable2)
#endif
#ifndef _INDEX_CONSERVED_QUANTITY_
#define _INDEX_CONSERVED_QUANTITY_(variable) (VARIABLE_1DLOOP,variable)
#endif
#ifndef _INDEX_VERTICAL_MOVEMENT_
#define _INDEX_VERTICAL_MOVEMENT_(variable) (VARIABLE_1DLOOP,variable)
#endif
#define _INDEX_EXTINCTION_ (VARIABLE_1DLOOP)

#else

! Not vectorized: RMBM subroutines operate one the local state only.

! Dummy argument and argument declaration for location specification.
#define LOCATION_ND LOCATION
#define DECLARE_LOCATION_ARG_ND LOCATION_TYPE,intent(in) :: LOCATION

! Beginning and end of spatial loop
#define _RMBM_LOOP_BEGIN_
#define _RMBM_LOOP_END_

! Dimensionality of generic space-dependent arguments.
#define ATTR_DIMENSIONS_0
#define ATTR_DIMENSIONS_1 ,dimension(:)
#define ATTR_DIMENSIONS_2 ,dimension(:,:)

! Expressions for indexing space-dependent RMBM variables defined on the full spatial domain.
#define _INDEX_ODE_(variable) (variable)
#define _INDEX_PPDD_(variable1,variable2) (variable1,variable2)
#define _INDEX_EXTINCTION_
#define _INDEX_CONSERVED_QUANTITY_(variable) (variable)
#define _INDEX_VERTICAL_MOVEMENT_(variable) (variable)

#endif

! For RMBM: standard arguments used in calling biogeochemical routines.
#define RMBM_ARGS_ND_IN root%environment,LOCATION_ND

! For BGC models: RMBM arguments to routines implemented by biogeochemical models.
#define RMBM_ARGS_ND environment,LOCATION_ND
#define RMBM_ARGS_DO_RHS RMBM_ARGS_ND,rhs
#define RMBM_ARGS_DO_PPDD RMBM_ARGS_ND,pp,dd
#define RMBM_ARGS_GET_EXTINCTION RMBM_ARGS_ND,extinction
#define RMBM_ARGS_GET_CONSERVED_QUANTITIES RMBM_ARGS_ND,sums
#define RMBM_ARGS_GET_SURFACE_EXCHANGE RMBM_ARGS_0D,flux

! For BGC models: Declaration of RMBM arguments to routines implemented by biogeochemical models.
#define DECLARE_RMBM_ARGS_ND type (type_environment),intent(inout) :: environment;DECLARE_LOCATION_ARG_ND
#define DECLARE_RMBM_ARGS_DO_RHS  DECLARE_RMBM_ARGS_ND;REALTYPE ATTR_DIMENSIONS_1,intent(inout) :: rhs
#define DECLARE_RMBM_ARGS_DO_PPDD DECLARE_RMBM_ARGS_ND;REALTYPE ATTR_DIMENSIONS_2,intent(inout) :: pp,dd
#define DECLARE_RMBM_ARGS_GET_EXTINCTION DECLARE_RMBM_ARGS_ND;REALTYPE ATTR_DIMENSIONS_0,intent(inout) :: extinction
#define DECLARE_RMBM_ARGS_GET_CONSERVED_QUANTITIES DECLARE_RMBM_ARGS_ND;REALTYPE ATTR_DIMENSIONS_1,intent(inout) :: sums
#define DECLARE_RMBM_ARGS_GET_SURFACE_EXCHANGE DECLARE_RMBM_ARGS_0D;REALTYPE,dimension(:),intent(inout) :: flux

! For BGC models: Expressions for setting space-dependent RMBM variables defined on the full spatial domain.
#define _SET_ODE_(variable,value) rhs _INDEX_ODE_(variable) = rhs _INDEX_ODE_(variable) + value
#define _SET_DD_(variable1,variable2,value) dd _INDEX_PPDD_(variable1,variable2) = dd _INDEX_PPDD_(variable1,variable2) + value
#define _SET_PP_(variable1,variable2,value) pp _INDEX_PPDD_(variable1,variable2) = pp _INDEX_PPDD_(variable1,variable2) + value
#define _SET_EXTINCTION_(value) extinction _INDEX_EXTINCTION_ = extinction _INDEX_EXTINCTION_ + value
#define _SET_CONSERVED_QUANTITY_(variable,value) sums _INDEX_CONSERVED_QUANTITY_(variable) = sums _INDEX_CONSERVED_QUANTITY_(variable) + value
#define _SET_VERTICAL_MOVEMENT_(variable,value) vertical_movement _INDEX_VERTICAL_MOVEMENT_(variable) = value
#define _SET_SURFACE_EXCHANGE_(variable,value) flux _INDEX_SURFACE_EXCHANGE_(variable) = value

! For BGC models: quick expressions for setting a single element in both the destriction and production matrix.
#define _SET_DD_SYM_(variable1,variable2,value) _SET_DD_(variable1,variable2,value);_SET_PP_(variable1,variable2,value)
#define _SET_PP_SYM_(variable1,variable2,value) _SET_PP_(variable1,variable2,value);_SET_DD_(variable1,variable2,value)

! For BGC models: read-only access to values of external dependencies
#define _GET_VAR_(variable) environment%var(variable)%data INDEX_LOCATION
#define _GET_VAR_HZ_(variable) environment%var_hz(variable)%data INDEX_LOCATION_HZ

! For BGC models: read-only access to state variable values
#ifdef RMBM_SINGLE_STATE_VARIABLE_ARRAY
#define _GET_STATE_(variable) environment%state _INDEX_STATE_(variable)
#else
#define _GET_STATE_(variable) environment%state(variable)%data INDEX_LOCATION
#endif

! For BGC models: write access to diagnostic variables
#ifdef RMBM_MANAGE_DIAGNOSTICS
#define _SET_DIAG_(index,value) environment%diag(index,LOCATION) = value
#define _SET_DIAG_HZ_(index,value) environment%diag_hz(index) = value
#else
#define _SET_DIAG_(index,value) environment%var(index)%data INDEX_LOCATION = value
#define _SET_DIAG_HZ_(index,value) environment%var_hz(index)%data INDEX_LOCATION_HZ = value
#endif

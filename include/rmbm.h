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
#define DEFINE_LOCATION_1DLOOP LOCATION_TYPE,intent(in) :: LOCATION_1DLOOP
#define ARG_LOCATION_1DLOOP ,LOCATION_1DLOOP
#else
#define DEFINE_LOCATION_1DLOOP
#define ARG_LOCATION_1DLOOP
#define LOCATION_1DLOOP
#define VARIABLE_1DLOOP LOCATION
#endif

! Define dimension attribute and index specifyer for horizontal (2D) fields.
#ifdef RMBM_HORIZONTAL_IS_SCALAR
#define INDEX_LOCATION_HZ
#define ATTR_LOCATION_DIMENSIONS_HZ
#else
#define INDEX_LOCATION_HZ (LOCATION_HZ)
#define ATTR_LOCATION_DIMENSIONS_HZ ,dimension(LOCATION_DIMENSIONS_HZ)
#endif

! Define dimension attribute and index specifyer for full 3D fields.
#define ATTR_LOCATION_DIMENSIONS ,dimension(LOCATION_DIMENSIONS)
#define INDEX_LOCATION (LOCATION)

#define _GET_STATE_(index) environment%state3d(index)%data INDEX_LOCATION
#define _GET_VAR_(index) environment%var3d(index)%data INDEX_LOCATION
#define _GET_VAR_HZ_(index) environment%var2d(index)%data INDEX_LOCATION_HZ
#define _SET_VAR_(index,value) environment%var3d(index)%data INDEX_LOCATION = value
#define _SET_VAR_HZ_(index,value) environment%var2d(index)%data INDEX_LOCATION_HZ = value

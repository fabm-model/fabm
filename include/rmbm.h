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
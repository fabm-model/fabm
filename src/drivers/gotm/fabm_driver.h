#define _FABM_DIMENSION_COUNT_ 1
#define _FABM_HORIZONTAL_IS_SCALAR_

! Variable name and dimension specifyer for full bio fields
#define LOCATION kk
#define LOCATION_DIMENSIONS :

! Let FABM manage (declare and allocate) arrays for diagnostic variables - slight performance gain.
#define _FABM_MANAGE_DIAGNOSTICS_

! Tell FABM that we want a 1D vectorized version
#define _FABM_USE_1D_LOOP_

! GOTM uses a different order of dimensions for rhs,pp,dd than expected by FABM.
! (FABM expects the vectorized spatial dimension to come first, which would benefit performance)
! Therefore, prescribe how to access spatially-localized rhs,pp,dd elements here.
! Note that this is only relevant if _FABM_USE_1D_LOOP_ is defined.
#define _INDEX_ODE_(i) (i,LOCATION-fabm_loop_start+1)
#define _INDEX_PPDD_(i,j) (i,j,LOCATION-fabm_loop_start+1)

! Include FABM preprocessor definitions.
! This *must* be done after the host-specific variables are defined (above),
! because these are used in fabm.h.
#include "fabm.h"

! Not used by FABM itself; only by GOTM-FABM driver gotm_fabm.F90
#define LOCATION_RANGE 0:kk


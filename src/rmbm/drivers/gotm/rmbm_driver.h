#define _RMBM_DIMENSION_COUNT_ 1
#define _RMBM_HORIZONTAL_IS_SCALAR_

! Variable name and dimension specifyer for full bio fields
#define LOCATION kk
#define LOCATION_DIMENSIONS :

! Let RMBM manage (declare and allocate) arrays for diagnostic variables - slight performance gain.
#define _RMBM_MANAGE_DIAGNOSTICS_

! Tell RMBM that we want a 1D vectorized version
#define _RMBM_USE_1D_LOOP_

! GOTM uses a different order of dimensions for rhs,pp,dd than expected by RMBM.
! (RMBM expects the vectorized spatial dimension to come first, which would benefit performance)
! Therefore, prescribe how to access spatially-localized rhs,pp,dd elements here.
! Note that this is only relevant if _RMBM_USE_1D_LOOP_ is defined.
#define _INDEX_ODE_(i) (i,LOCATION-rmbm_loop_start+1)
#define _INDEX_PPDD_(i,j) (i,j,LOCATION-rmbm_loop_start+1)

! Include RMBM preprocessor definitions.
! This *must* be done after the host-specific variables are defined (above),
! because these are used in rmbm.h.
#include "rmbm.h"

! Not used by RMBM itself; only by GOTM-RMBM driver gotm_rmbm.F90
#define LOCATION_RANGE 0:kk


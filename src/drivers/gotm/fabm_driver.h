#define _FABM_DIMENSION_COUNT_ 1
#define _FABM_DEPTH_DIMENSION_INDEX_ 1
#define _FABM_VECTORIZED_DIMENSION_INDEX_ 1

! GOTM uses a different order of dimensions for rhs,pp,dd than expected by FABM.
! (FABM expects the vectorized spatial dimension to come first, which would benefit performance)
! Therefore, prescribe how to access spatially-localized rhs,pp,dd elements here.
! Note that this is only relevant if _FABM_VECTORIZED_DIMENSION_INDEX_ is defined.
#define _INDEX_ODE_(i) (i,i__-fabm_loop_start+1)
#define _INDEX_ODE_1D_(i) (i,1:fabm_loop_stop-fabm_loop_start+1)
#define _INDEX_PPDD_(i,j) (i,j,i__-fabm_loop_start+1)

! Include FABM preprocessor definitions.
! This *must* be done after the host-specific variables are defined (above),
! because these are used in fabm.h.
#include "fabm.h"

! Not used by FABM itself; only by GOTM-FABM driver gotm_fabm.F90
#define _LOCATION_RANGE_ 0:i__


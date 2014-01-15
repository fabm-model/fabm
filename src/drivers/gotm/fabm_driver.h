#define _FABM_DIMENSION_COUNT_ 1
#define _FABM_DEPTH_DIMENSION_INDEX_ 1
#define _FABM_VECTORIZED_DIMENSION_INDEX_ 1

#define _FABM_VERTICAL_BOTTOM_TO_SURFACE_

! Include FABM preprocessor definitions.
! This *must* be done after the host-specific variables are defined (above),
! because these are used in fabm.h.
#include "fabm.h"

! Not used by FABM itself; only by GOTM-FABM driver gotm_fabm.F90
#define _LOCATION_RANGE_ 0:i__


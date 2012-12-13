#define _FABM_DIMENSION_COUNT_ 1
#define _FABM_HORIZONTAL_IS_SCALAR_

! Variable name and dimension specifier for full bio fields
#define _LOCATION_ kk
#define _LOCATION_DIMENSIONS_ :

! Let FABM manage (declare and allocate) arrays for diagnostic variables
#define _FABM_MANAGE_DIAGNOSTICS_

! Include FABM preprocessor definitions.
! This *must* be done after the host-specific variables are defined (above),
! because these are used in fabm.h.
#include "fabm.h"

! Not used by FABM itself; only by GLM-FABM driver glm_fabm.F90
#define _LOCATION_RANGE_ 1:kk


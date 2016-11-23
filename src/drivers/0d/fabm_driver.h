#undef REALTYPE
#undef _ZERO_
#undef _ONE_
#define REALTYPE real(rk)
#define _ZERO_ 0._rk
#define _ONE_  1._rk

#define _FABM_DIMENSION_COUNT_ 0

! Include FABM preprocessor definitions.
! This *must* be done after the host-specific variables are defined (above),
! because these are used in fabm.h.
#include "fabm.h"


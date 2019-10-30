#define _FABM_DIMENSION_COUNT_ 3
#define _FABM_DEPTH_DIMENSION_INDEX_ 3
#define _FABM_VECTORIZED_DIMENSION_INDEX_ 1

#define _FABM_MASK_TYPE_ real(rke)
#define _FABM_UNMASKED_VALUE_ 1

! Specify that the vertical index of the bottom cell is variable (kmax depends on i,j)
#define _FABM_BOTTOM_INDEX_ -1

#include "fabm.h"


#ifndef _FABM_DIMENSION_COUNT_
#  error "_FABM_DIMENSION_COUNT_ must be defined."
#endif

#if _FABM_DIMENSION_COUNT_==1
#  define _FABM_VECTORIZED_DIMENSION_INDEX_ 1
#  define _FABM_CONTIGUOUS_
#endif

#include "fabm.h"


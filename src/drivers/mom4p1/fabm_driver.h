#define _FABM_DIMENSION_COUNT_ 3

! Variable name and dimension specifyer for full bio fields
#define _LOCATION_ ii,jj,kk
#define _LOCATION_DIMENSIONS_ :,:,:

! Variable name and dimension specifyer for horizontal-only bio fields
#define _LOCATION_HZ_ ii,jj
#define _LOCATION_DIMENSIONS_HZ_ :,:

! Divide location variables over a single iterator variable and remaining variables.
#define _VARIABLE_1DLOOP_ ii
#define _LOCATION_1DLOOP_ jj,kk

!#define _FABM_USE_1D_LOOP_

#include "fabm.h"


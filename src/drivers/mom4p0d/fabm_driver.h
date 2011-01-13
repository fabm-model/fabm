#define RMBM_DIMENSIONS 3

! Variable name and dimension specifyer for full bio fields
#define LOCATION ii,jj,kk
#define LOCATION_DIMENSIONS :,:,:

! Variable name and dimension specifyer for horizontal-only bio fields
#define LOCATION_HZ ii,jj
#define LOCATION_DIMENSIONS_HZ :,:

! Divide location variables over a single iterator variable and remaining variables.
#define VARIABLE_1DLOOP ii
#define LOCATION_1DLOOP jj,kk

#include "rmbm.h"


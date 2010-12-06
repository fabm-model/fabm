#define RMBM_DIMENSIONS 1
#define RMBM_HORIZONTAL_IS_SCALAR

! Variable name and dimension specifyer for full bio fields
#define LOCATION kk
#define LOCATION_DIMENSIONS :

! Variable name and dimension specifyer for horizontal-only bio fields
#define LOCATION_HZ
#define LOCATION_DIMENSIONS_HZ

! Include RMBM preprocessor definitions.
! This *must* be done after the host-specific variables are defined (above),
! because these are used in rmbm.h.
#include "rmbm.h"

! Not used by RMBM itself; only by GOTM-RMBM driver gotm_rmbm.F90
#define LOCATION_RANGE 0:kk

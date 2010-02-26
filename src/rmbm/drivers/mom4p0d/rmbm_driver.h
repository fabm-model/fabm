#include "rmbm.h"

#define LOCATIONTYPE integer
#define LOCATION ii,jj,kk
#define LOCATIONDIMENSIONS :,:,:
#define LOCATION2D ii,jj
#define LOCATION2DDIMENSIONS :,:

#define LOCATION_1DLOOP jj,kk
#define VARIABLE_1DLOOP ii
#define LENGTH_1DLOOP iec-isc+1

#define LOOP1D_USE use ocean_tpm_util_mod, only: isc, iec

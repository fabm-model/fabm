#include "rmbm.h"

#define LOCATIONTYPE integer
#define LOCATION ii,jj,kk
#define LOCATIONDIMENSIONS :,:,:
#define LOCATION2D ii,jj
#define LOCATION2DDIMENSIONS :,:

#define LOCATION_1DLOOP jj,kk
#define SLICE1D ii,:
#define LOOP1D ii=1,iec-isc+1
#define VARIABLE_1DLOOP ii

#define LOOP1D_USE use ocean_tpm_util_mod, only: isc, iec

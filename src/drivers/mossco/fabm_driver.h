#define _FABM_DIMENSION_COUNT_ 3
#define _FABM_DEPTH_DIMENSION_INDEX_ 3

! Include FABM preprocessor definitions.
! This *must* be done after the host-specific variables are defined (above),
! because these are used in fabm.h.
#include "fabm.h"

! Taken from the General Ocean Turbulence Model:

#define PATH_MAX	255

#define stderr		0
#define stdout		6

#define STDOUT write(stdout,*)
#define STDERR write(stderr,*)
#define LEVEL0 STDERR
#define LEVEL1 STDERR '   ',
#define LEVEL2 STDERR '       ',
#define LEVEL3 STDERR '           ',
#define LEVEL4 STDERR '               ',
#define FATAL  STDERR 'FATAL ERROR: ',

#define LINE "------------------------------------------------------------------------"

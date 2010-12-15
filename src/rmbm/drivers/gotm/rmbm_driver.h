#define RMBM_DIMENSIONS 1
#define RMBM_HORIZONTAL_IS_SCALAR

! Variable name and dimension specifyer for full bio fields
#define LOCATION kk
#define LOCATION_DIMENSIONS :

! Let RMBM manage (declare and allocate) arrays for diagnsotic variables - slight performance gain.
#define RMBM_MANAGE_DIAGNOSTICS

! Tell RMBM that we will store all state variables in a single array - slight performance gain.
! Also tell RMBM how to retrieve a spatially-localized state variable values from that array.
#define RMBM_SINGLE_STATE_VARIABLE_ARRAY
#define _INDEX_STATE_(variable) (variable,LOCATION)

! Tell RMBM that we want a 1D vectorized version
#define RMBM_USE_1D_LOOP

! GOTM uses a different order of dimensions for rhs,pp,dd than expected by RMBM.
! (it expects the vectorized spatial dimension to come first, which would benefit performance)
! Therefore, prescribe how to access spatially-localized rhs,pp,dd elements here.
! Note that this is only relevant if RMBM_USE_1D_LOOP is defined.
#define _INDEX_ODE_(i) (i,LOCATION)
#define _INDEX_PPDD_(i,j) (i,j,LOCATION)

! Include RMBM preprocessor definitions.
! This *must* be done after the host-specific variables are defined (above),
! because these are used in rmbm.h.
#include "rmbm.h"

! Not used by RMBM itself; only by GOTM-RMBM driver gotm_rmbm.F90
#define LOCATION_RANGE 0:kk


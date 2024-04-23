! SPDX-FileCopyrightText: 2021-2024 Helmholtz-Zentrum hereon GmbH
! SPDX-FileCopyrightText: 2017-2021 Helmholtz-Zentrum Geesthacht GmbH
! SPDX-License-Identifier: CC0-1.0
! SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
! SPDX-FileContributor: Richard Hofmeister

#define _FABM_DIMENSION_COUNT_ 2
#define _FABM_DEPTH_DIMENSION_INDEX_ 1
#define _FABM_VECTORIZED_DIMENSION_INDEX_ 1

#define _FABM_VERTICAL_BOTTOM_TO_SURFACE_
#define _FABM_BOTTOM_INDEX_ -1

#define _FABM_MASK_TYPE_ integer
#define _FABM_MASKED_VALUE_ 1

! For performance reasons it would be advisable to define _FABM_HORIZONTAL_MASK_
! which is sufficient to mask dry elements.  It currently does break the 
! code however, which still employs a 2D mask to mask elements below the bottom index
! within a water column.  We need a discussion on whether it is more performant to 
! do the calculations below bottom index nonetheless ... 
   
#include "fabm.h"


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

#include "fabm.h"


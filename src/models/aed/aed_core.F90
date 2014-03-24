!###############################################################################
!#                                                                             #
!# aed_core.F90                                                                #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created June 2013                                                           #
!#                                                                             #
!###############################################################################

#include "aed.h"

MODULE aed_core
   USE fabm_types

   IMPLICIT NONE

! Re-export fabm bits
!

!CONSTANTS
   AED_REAL,PARAMETER :: zero_ = 0., one_ = 1.
   AED_REAL,PARAMETER :: secs_per_day = 86400.
   AED_REAL,PARAMETER :: misval_ = -9999.

!===============================================================================
!CONTAINS


END MODULE aed_core

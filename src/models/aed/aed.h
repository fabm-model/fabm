!###############################################################################
!#                                                                             #
!# aed.h                                                                       #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!###############################################################################
#ifndef _AED_H_
#define _AED_H_

#define AED_VERSION  '0.35'

#define MAX_MODELS 40

!# aed_phytoplankton constants
#define MAX_PHYTO_TYPES 256
#define MAX_ZOOP_TYPES  256
#define MAX_ZOOP_PREY    10
#define MAX_PATHO_TYPES 256

!# for aed_geochemistry
#define MAX_GC_COMPONENTS 20
#define MAX_GC_MINERALS   20

#define MISVAL -9999.

#define INP_LINE_LEN 512
#define STR_LEN       32
#define REACTION_START_CH  '['
#define REACTION_END_CH    ']'
#define PLUS               '+'
#define MINUS              '-'
#define EQUALS             '='
#define L_PAREN            '('
#define R_PAREN            ')'

#include "fabm_driver.h"
#define AED_REAL real(rk)
#define namlst configunit

#define DOUBLETYPE double precision

#endif

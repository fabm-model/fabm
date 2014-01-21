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

#define AED_VERSION  '0.32'

#include "fabm_driver.h"

! aed_phytoplankton constants
#define MAX_PHYTO_TYPES 256
#define MAX_ZOOP_TYPES  256
#define MAX_ZOOP_PREY    10
#define MAX_PATHO_TYPES 256

! for aed_geochemistry
#define MAX_GC_COMPONENTS 20
#define MAX_GC_MINERALS   20

#define MISVAL -9999.

!#ifdef SINGLE
!#define REALTYPE real
!#define r_wq 4
!#else
!#define REALTYPE double precision
!#define IFIX IDINT
!#define AMOD DMOD
!#define ALOG10 DLOG10
!#define r_wq 8
!#endif
#define DOUBLETYPE double precision

!#define stdin  5
!#define stdout 6
!#define stderr 0

#define INP_LINE_LEN 512
#define STR_LEN       32
#define REACTION_START_CH  '['
#define REACTION_END_CH    ']'
#define PLUS               '+'
#define MINUS              '-'
#define EQUALS             '='
#define L_PAREN            '('
#define R_PAREN            ')'

#define DO_DEBUG 1
#if DO_DEBUG
# define ASSERT(expr)  if (.not.(expr)) print *,"assert failed in ",__FILE__," at line ",__LINE__
# define _CheckAllocStatus(status) if (status /= 0) print *,"(de)allocation error in ",__FILE__," at line ",__LINE__
# define CheckAllocStatus(status,a,b) if (status /= 0) print *,"(de)allocation error in ",__FILE__," at line ",__LINE__
# define _CheckFileIOStatus(status) if (status /= 0) print *,"file io error in ",__FILE__," at line ",__LINE__
# define CheckFileIOStatus(a,b,status,c,d) if (status /= 0) print *,"file io error in ",__FILE__," at line ",__LINE__
#else
# define ASSERT(expr)
# define _CheckAllocStatus(status)
# define CheckAllocStatus(status,a,b)
# define _CheckFileIOStatus(status)
# define CheckFileIOStatus(a,b,status,c,d)
#endif

#endif

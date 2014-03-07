!###############################################################################
!#                                                                             #
!# aed_totals.F90                                                              #
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
!# Created July 2012                                                           #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_totals
!-------------------------------------------------------------------------------
! aed_totals --- totals biogeochemical model
!
! The AED module totals contains only diagnostic variables to provide
! totals of other variables (eg tss)
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_type_totals
!
   TYPE,extends(type_base_model) :: aed_type_totals
!     Variable identifiers
      TYPE (type_diagnostic_variable_id) :: id_totals_tn, id_totals_tp, id_totals_toc,  &
                                            id_totals_tss, id_totals_turbidity
      TYPE (type_state_variable_id),ALLOCATABLE :: id_dep_tn(:), id_dep_tp(:), id_dep_toc(:),  &
                                              id_dep_tss(:)
      AED_REAL,ALLOCATABLE               :: turbidity(:)


!     Model parameters
      AED_REAL :: Fsed_dic,Ksed_dic,theta_sed_dic
      LOGICAL  :: use_oxy,use_dic

      CONTAINS      ! Model Methods
        PROCEDURE :: initialize => aed_init_totals
        PROCEDURE :: do         => aed_totals_do
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_totals(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_type_totals),TARGET,INTENT(inout) :: self

!
!LOCALS
   INTEGER :: status

   INTEGER           :: i, num_tn,num_tp,num_toc,num_tss
   CHARACTER(len=40) :: tn(100), tp(100), toc(100), tss(100)
   AED_REAL          :: turbidity(100)

   NAMELIST /aed_totals/ tn,tp,toc,tss,turbidity
!
!-------------------------------------------------------------------------------
!BEGIN
   tn = '' ; tp = '' ; toc = '' ; tss = '' ; turbidity = MISVAL

   ! Read the namelist
   read(namlst,nml=aed_totals,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_totals'

   DO i=1,100 ; IF (tn(i)  .EQ. '' ) THEN ; num_tn  = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tp(i)  .EQ. '' ) THEN ; num_tp  = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (toc(i) .EQ. '' ) THEN ; num_toc = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (tss(i) .EQ. '' ) THEN ; num_tss = i-1 ; EXIT ; ENDIF ; ENDDO

!  print *,"totl tn = ",num_tn," num_tp = ",num_tp," num toc = ",num_toc," num_tss = ",num_tss
   ALLOCATE(self%id_dep_tn(num_tn))
   ALLOCATE(self%id_dep_tp(num_tp))
   ALLOCATE(self%id_dep_toc(num_toc))
   ALLOCATE(self%id_dep_tss(num_tss))
   ALLOCATE(self%turbidity(num_tss))

   ! Register external state variable dependencies
   DO i=1,num_tn  ; call self%register_state_dependency(self%id_dep_tn(i),tn(i))  ; ENDDO
   DO i=1,num_tp  ; call self%register_state_dependency(self%id_dep_tp(i),tp(i))  ; ENDDO
   DO i=1,num_toc ; call self%register_state_dependency(self%id_dep_toc(i),toc(i)) ; ENDDO
   DO i=1,num_tss ; call self%register_state_dependency(self%id_dep_tss(i),tss(i)) ; ENDDO
   self%turbidity = turbidity(1:num_tss)

   ! Register diagnostic variables
   CALL self%register_diagnostic_variable(self%id_totals_tn,'tn',               &
                     'mmol/m**2/d', 'Filterable reactive totals')

   CALL self%register_diagnostic_variable(self%id_totals_tp,'tp',               &
                     'mmol/m**2/d', 'Filterable reactive totals')

   CALL self%register_diagnostic_variable(self%id_totals_toc,'toc',             &
                     'mmol/m**2/d', 'Filterable reactive totals')

   CALL self%register_diagnostic_variable(self%id_totals_tss,'tss',             &
                     'mmol/m**2/d', 'Filterable reactive totals')

   CALL self%register_diagnostic_variable(self%id_totals_turbidity,'turbidity', &
                     'mmol/m**2/d', 'Filterable reactive totals')
END SUBROUTINE aed_init_totals
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_totals_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_totals model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_totals),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
   INTEGER :: i,count
   AED_REAL :: val, tot, tot2

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   tot = 0.
   count = ubound(self%id_dep_tn,1)
   DO i=1,count ; _GET_(self%id_dep_tn(i),val) ; tot = tot + val ; ENDDO
   _SET_DIAGNOSTIC_(self%id_totals_tn, tot)

   tot = 0.
   count = ubound(self%id_dep_tp,1)
   DO i=1,count ; _GET_(self%id_dep_tp(i),val) ; tot = tot + val ; ENDDO
   _SET_DIAGNOSTIC_(self%id_totals_tp, tot)

   tot = 0.
   count = ubound(self%id_dep_toc,1)
   DO i=1,count ; _GET_(self%id_dep_toc(i),val) ; tot = tot + val ; ENDDO
   _SET_DIAGNOSTIC_(self%id_totals_toc, tot)

   tot = 0.
   tot2 = 0.
   count = ubound(self%id_dep_tss,1)
   DO i=1,count
      _GET_(self%id_dep_tss(i),val)
      tot = tot + val
      tot2 = tot2 + val*(self%turbidity(i))
   ENDDO

   _SET_DIAGNOSTIC_(self%id_totals_tss, tot)
   _SET_DIAGNOSTIC_(self%id_totals_turbidity, tot2)

   ! Leave spatial loops (if any)
   _LOOP_END_
END SUBROUTINE aed_totals_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_totals

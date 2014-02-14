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
   USE fabm_types

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC type_aed_totals
!
   TYPE,extends(type_base_model) :: type_aed_totals
!     Variable identifiers
      type (type_diagnostic_variable_id)        :: id_totals_tn, id_totals_tp, id_totals_toc,  &
                                              id_totals_tss, id_totals_turbidity
      type (type_state_variable_id),ALLOCATABLE :: id_dep_tn(:), id_dep_tp(:), id_dep_toc(:),  &
                                              id_dep_tss(:)
      real(rk),ALLOCATABLE                 :: turbidity(:)


!     Model parameters
      real(rk) :: Fsed_dic,Ksed_dic,theta_sed_dic
      LOGICAL  :: use_oxy,use_dic

      CONTAINS      ! Model Methods
        procedure :: initialize               => aed_totals_init
        procedure :: get_conserved_quantities => aed_totals_get_conserved_quantities
   END TYPE


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_totals_init(self,configunit)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (type_aed_totals),TARGET,INTENT(INOUT) :: self
   INTEGER,INTENT(in)                           :: configunit
!
!LOCALS

   INTEGER           :: i, num_tn,num_tp,num_toc,num_tss
   CHARACTER(len=40) :: tn(100), tp(100), toc(100), tss(100)
   real(rk)          :: turbidity(100)

   NAMELIST /aed_totals/ tn,tp,toc,tss,turbidity
!
!-------------------------------------------------------------------------------
!BEGIN
   tn = '' ; tp = '' ; toc = '' ; tss = '' ; turbidity = MISVAL

   ! Read the namelist
   read(configunit,nml=aed_totals,err=99)

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
   call self%register_diagnostic_variable(self%id_totals_tn,'aed_totals_tn',               &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated)

   call self%register_diagnostic_variable(self%id_totals_tp,'aed_totals_tp',               &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated)

   call self%register_diagnostic_variable(self%id_totals_toc,'aed_totals_toc',             &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated)

   call self%register_diagnostic_variable(self%id_totals_tss,'aed_totals_tss',             &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated)

   call self%register_diagnostic_variable(self%id_totals_turbidity,'aed_totals_turbidity', &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated)


   RETURN

99 CALL self%fatal_error('aed_totals_init','Error reading namelist aed_totals')

END SUBROUTINE aed_totals_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_totals_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_totals model
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_totals),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   INTEGER :: i,count
   real(rk) :: val, tot, tot2

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

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
   _FABM_LOOP_END_
END SUBROUTINE aed_totals_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_totals

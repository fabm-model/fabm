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

#ifdef _FABM_F2003_

#include "fabm_driver.h"
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
   USE fabm_driver

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC type_aed_totals, aed_totals_create
!
   TYPE,extends(type_base_model) :: type_aed_totals
!     Variable identifiers
      _TYPE_DIAGNOSTIC_VARIABLE_ID_        :: id_totals_tn, id_totals_tp, id_totals_toc,  &
                                              id_totals_tss, id_totals_turbidity
      _TYPE_STATE_VARIABLE_ID_,ALLOCATABLE :: id_dep_tn(:), id_dep_tp(:), id_dep_toc(:),  &
                                              id_dep_tss(:)
      REALTYPE,ALLOCATABLE                 :: turbidity(:)


!     Model parameters
      REALTYPE :: Fsed_dic,Ksed_dic,theta_sed_dic
      LOGICAL  :: use_oxy,use_dic

      CONTAINS      ! Model Methods
!       procedure :: initialize               => aed_totals_init
        procedure :: get_conserved_quantities => aed_totals_get_conserved_quantities
   END TYPE


!===============================================================================
CONTAINS


!###############################################################################
FUNCTION aed_totals_create(namlst,name,parent) RESULT(self)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CHARACTER(len=*),INTENT(in) :: name
   _CLASS_ (type_model_info),TARGET,INTENT(inout) :: parent
!
!LOCALS
   _CLASS_ (type_aed_totals),POINTER :: self

   INTEGER           :: i, num_tn,num_tp,num_toc,num_tss
   CHARACTER(len=40) :: tn(100), tp(100), toc(100), tss(100)
   REALTYPE          :: turbidity(100)

   NAMELIST /aed_totals/ tn,tp,toc,tss,turbidity
!
!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(self)
   CALL self%initialize(name,parent)

   tn = '' ; tp = '' ; toc = '' ; tss = '' ; turbidity = MISVAL

   ! Read the namelist
   read(namlst,nml=aed_totals,err=99)

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
   DO i=1,num_tn  ; self%id_dep_tn(i)  = self%register_state_dependency(tn(i))  ; ENDDO
   DO i=1,num_tp  ; self%id_dep_tp(i)  = self%register_state_dependency(tp(i))  ; ENDDO
   DO i=1,num_toc ; self%id_dep_toc(i) = self%register_state_dependency(toc(i)) ; ENDDO
   DO i=1,num_tss ; self%id_dep_tss(i) = self%register_state_dependency(tss(i)) ; ENDDO
   self%turbidity = turbidity(1:num_tss)

   ! Register diagnostic variables
   self%id_totals_tn = self%register_diagnostic_variable('aed_totals_tn',               &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated, shape=shape_hz)

   self%id_totals_tp = self%register_diagnostic_variable('aed_totals_tp',               &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated, shape=shape_hz)

   self%id_totals_toc = self%register_diagnostic_variable('aed_totals_toc',             &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated, shape=shape_hz)

   self%id_totals_tss = self%register_diagnostic_variable('aed_totals_tss',             &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated, shape=shape_hz)

   self%id_totals_turbidity = self%register_diagnostic_variable('aed_totals_turbidity', &
                     'mmol/m**2/d', 'Filterable reactive totals',                       &
                     time_treatment=time_treatment_step_integrated, shape=shape_hz)


   RETURN

99 CALL fatal_error('aed_totals_init','Error reading namelist aed_totals')

END FUNCTION aed_totals_create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_totals_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_totals model
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_totals),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   INTEGER :: i,count
   REALTYPE :: val, tot, tot2

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   tot = 0.
   count = ubound(self%id_dep_tn,1)
   DO i=1,count ; _GET_STATE_(self%id_dep_tn(i),val) ; tot = tot + val ; ENDDO
   _SET_DIAG_(self%id_totals_tn, tot)

   tot = 0.
   count = ubound(self%id_dep_tp,1)
   DO i=1,count ; _GET_STATE_(self%id_dep_tp(i),val) ; tot = tot + val ; ENDDO
   _SET_DIAG_(self%id_totals_tp, tot)

   tot = 0.
   count = ubound(self%id_dep_toc,1)
   DO i=1,count ; _GET_STATE_(self%id_dep_toc(i),val) ; tot = tot + val ; ENDDO
   _SET_DIAG_(self%id_totals_toc, tot)

   tot = 0.
   tot2 = 0.
   count = ubound(self%id_dep_tss,1)
   DO i=1,count
      _GET_STATE_(self%id_dep_tss(i),val)
      tot = tot + val
      tot2 = tot2 + val*(self%turbidity(i))
   ENDDO

   _SET_DIAG_(self%id_totals_tss, tot)
   _SET_DIAG_(self%id_totals_turbidity, tot2)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_
END SUBROUTINE aed_totals_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_totals
#endif
